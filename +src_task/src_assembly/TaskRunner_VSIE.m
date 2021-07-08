classdef TaskRunner_VSIE < TaskRunner_Base
    
    
    % ------------------------------------------------------------- %

    properties
        
        % get coil and body models
        coil
        scatterer
         
        % assemble VSIE system
        assembly_vsie
        
        % iterative solver 
        solver
        
        % temporary var. for final currents
        Jcb_sol
        
        sol_time 
        rel_res
        res_vec = {}
        
        % port network parameters
        NP     
    end
    
    % ------------------------------------------------------------- %
    
    methods
        
        % define Task_VSIE constructor
        function [obj] = TaskRunner_VSIE(task_settings)
            
            % call base class constructor
            obj = obj@TaskRunner_Base(task_settings);
            
            % set coil loader and load a coil
            obj.coil_loader = src_geo.Coil_Loader(task_settings.sie);
            
            % load coil and body models
            obj.load_coil_model_();
            obj.load_scatterer_();
            
            % set up task assembler
            obj.assembly_vsie = src_assembly.Assembly_VSIE(task_settings, obj.coil,...
                                                           obj.scatterer, obj.dims);
        end
        
        % ------------------------------------------------------------- %

        % solve the problem
        function run(obj)
            
            % get working frequencies
            [freqs, N_freq] = obj.get_working_frequencies_();
            
            % run frequency sweep 
            for  i_freq = 1:N_freq
                
                % get current frequency
                freq = freqs(i_freq);
                
                % if the system is already assembled - do not reassemble
                if isempty(obj.solver) || (freq ~= freqs(1))
                    
                    % assemble VIE system
                    obj.assemble_VSIE_system_(freq);
                    
                    % setup iterative solver
                    obj.solver = obj.assembly_vsie.setup_solver();
                end
                
                if isempty(obj.Jcb_sol)
                    % initialize solution
                    obj.init_solution_();
                end
                
                obj.compute_currents_(i_freq);
                
                % compute electric, magnetic fields
                fprintf('Computing Electric field ... \n\n');

                obj.compute_E_total_(i_freq);
                
                fprintf('Computing Magnetic field ... \n\n');
                
                obj.compute_H_total_(i_freq);

            end  
            
            % update solution object
            obj.update_solution_();
            
            % compute absorbed power, b1+,b1-
            obj.compute_P_abs_();
            obj.compute_SAR_();
            obj.compute_b1p_();
            obj.compute_b1m_();
            
        end
        
    end
    
    % ------------------------------------------------------------- %

    methods (Access = private)
        
        % ------------------------------------------------------------- %
       
        function [] = load_coil_model_(obj)
            % private method loads coil to task runner
            obj.coil = obj.coil_loader.load_coil();
            
            % get coil dimensions
            obj.dims.get_coil_dims(obj.coil);
        end
        
        % ------------------------------------------------------------- %
        
        function [] = assemble_VSIE_system_(obj, freq)
            obj.assembly_vsie.assemble_system(freq);
        end
        
        % ------------------------------------------------------------- %
        
        function load_scatterer_(obj)
            % function load_scatterer(obj, freq) 
            
            % get current frequency
            [freqs, ~] = obj.get_working_frequencies_();
                        
            % load RHBM structure from file
            RHBM = load(obj.task_settings_.vie.BodyModel);
            
            % init scatterer object
            obj.scatterer = src_scatterer.Scatterer(RHBM, obj.dims, freqs(1));
            
        end

        % ------------------------------------------------------------- %
        
        function init_solution_(obj)
            % function initializes solution object
            
            % get dimensions of the problem
            dims = obj.dims;
            
            % initialize buffer
            obj.Jcb_sol = zeros(dims.N_sie+dims.N_scat*dims.ql, dims.N_feeds,dims.N_freqs);
            
            % initialize solution object for VSIE problem
            obj.task_solution_.init_VSIE_solution(dims);
        end
        
        % ------------------------------------------------------------- %
        
        function [] = compute_currents_(obj, freq_num)
            
            w_src_I  = obj.assembly_vsie.w_src_I;
            w_coil_V = obj.assembly_vsie.w_coil_V;
            Zw_coil = obj.assembly_vsie.Zw_coil;
            Zw_src  = obj.assembly_vsie.Zw_src;
            Yw_coil = obj.assembly_vsie.Yw_coil;
            Fcoil = obj.assembly_vsie.Fcoil_;
                        
            % loop over rhs vectors
                for feed_port = 1:obj.dims.N_feeds
                    
                    fprintf('Solving VSIE problem for Port %d \n', feed_port);

                    tic
                    
                    switch obj.task_settings_.vsie.Solver_mode
                        
                        case 'Decoupled'
                            rhs       = obj.assembly_vsie.rhs_b;
                            ini_guess = zeros(obj.dims.N_scat * obj.dims.ql,1);

                        case 'Coupled'
                            
                            rhs = obj.assembly_vsie.rhs;
                            
                            ini_Jc = obj.assembly_vsie.Zcoil_inv_ * rhs(1:obj.dims.N_sie, feed_port);
                            ini_guess = [ini_Jc; zeros(obj.dims.N_scat * obj.dims.ql,1)];
 
                    end
                    
                    % run solver 
                    [Jcb, relative_res, res_vector] = obj.solver.run(rhs, ini_guess, feed_port);

                    solve_toc = toc;

                    % store results in temp tensors
                    obj.Jcb_sol(:,feed_port, freq_num) = Jcb;
                    obj.sol_time(feed_port,freq_num)   = solve_toc;
                    obj.rel_res(:,feed_port, freq_num) = relative_res;
                    obj.res_vec(feed_port, freq_num)   = {res_vector};                   
                    
                    fprintf('VSIE problem solved in %d iterations\n', length(res_vector));
                    fprintf('Relative residual: %e \n', relative_res);
                    fprintf('Elapsed time: %4.2f seconds \n\n', solve_toc);

                end
                
                Jcoil = squeeze(obj.Jcb_sol(1:obj.dims.N_sie,:,freq_num));
                
                % the actual voltage applied to excitation ports (Coil) 
                Icoil = - Fcoil.' * Jcoil;
                Vcoil = obj.assembly_vsie.rhs_cp + obj.assembly_vsie.Z_L_hat * Fcoil.' * Jcoil;
                
                Isrc  = w_src_I * Icoil + Yw_coil * Vcoil;
                Vsrc  = w_coil_V * Vcoil + Zw_coil * Icoil + Zw_src * Isrc;
                
                % find NP seeing from the voltage sourse (i.e. coil + external circuit + ref. impedance)
                obj.NP.src.Y(:,:,freq_num) = Isrc / Vsrc;
                obj.NP.src.Z(:,:,freq_num) = y2z(obj.NP.src.Y(:,:,freq_num));
                obj.NP.src.S(:,:,freq_num) = y2s(obj.NP.src.Y(:,:,freq_num));
                
                % find NP seeing from the coil ports
                obj.NP.coil.Y(:,:,freq_num) = Icoil / Vcoil;
                obj.NP.coil.Z(:,:,freq_num) = y2z(obj.NP.coil.Y(:,:,freq_num));
                obj.NP.coil.S(:,:,freq_num) = y2s(obj.NP.coil.Y(:,:,freq_num));
           
        end
        
        % ------------------------------------------------------------- %
        
        function update_solution_(obj)
            
            % get dimensions of the problem
            dims = obj.dims;
            
            [freqs, N_freq] = obj.get_working_frequencies_();
            idx_ql   = obj.scatterer.index_vie.index_ql(dims.ql, dims.Nvox_vie);
            
            Jc        = obj.Jcb_sol(1:dims.N_sie,:,:);
            Jb_vec    = obj.Jcb_sol(dims.N_sie+1:end,:,:);
            Jb_tensor = zeros([dims.vie, dims.N_feeds, N_freq]);

            for i = 1:dims.N_feeds
                for j = 1:N_freq
                    
                    N_shift = dims.Nvox_vie * dims.ql;
                    idx = idx_ql + (i-1) * N_shift + ...
                          (j-1) * N_shift * dims.N_feeds;
                    Jb_tensor(idx) = Jb_vec(:,i,j);
                end
            end    
            
            % store solution of the VSIE problem in task_solution
            obj.task_solution_.set_surface_currents(Jc);
            obj.task_solution_.set_plrz_currents_vec(Jb_vec);
            obj.task_solution_.set_plrz_currents_tensor(Jb_tensor);
            
            % store relative residual and res. vector
            obj.task_solution_.set_relative_residual(obj.rel_res);
            obj.task_solution_.set_residual_vector(obj.res_vec);
            obj.task_solution_.set_lapse_time(obj.sol_time);
            
            % set network parameters
            obj.task_solution_.set_Network_param(obj.NP);
            
            % store total electric and magnetic fields
            obj.task_solution_.set_E_total(obj.Etot)
            obj.task_solution_.set_H_total(obj.Htot)
            
            % set the relative permittivity and conductivity
            obj.task_solution_.set_epsilon_r(obj.scatterer.prop_vie.epsilon_r);
            obj.task_solution_.set_sigma_e(obj.scatterer.prop_vie.sigma_e);
            
            % set the frequencies
            obj.task_solution_.set_freqs(freqs);
                        
        end
        
        % ------------------------------------------------------------- %
        
        function [] = compute_E_total_(obj, i_freq)
        % compute total electric fields 
        
            Jcb = obj.Jcb_sol(obj.dims.N_sie+1:end,:,:);
            
            % get tissue currents 
            obj.compute_E_total_base_(Jcb, obj.assembly_vsie.mvp,...
                                       i_freq)
        
        end
        

        % ------------------------------------------------------------- %
        
        function [] = compute_H_total_(obj, i_freq)
                        
            % computes the magnetic field
            coupling_mode = obj.task_settings_.vsie.Coupling_Mode;
            
            Jc  = obj.Jcb_sol(1:obj.dims.N_sie,:, i_freq);
            Jb  = obj.Jcb_sol(obj.dims.N_sie+1:end,:, i_freq);
            mvp = obj.assembly_vsie.mvp;
            
            switch coupling_mode
                case 'pFFT'
                    obj.compute_Htot_pfft_(Jc, Jb, mvp, i_freq);
                    
                case 'Basis'
                    obj.compute_Htot_basis_(Jc, Jb, mvp, i_freq);
                    
                case 'Dense'
                    obj.compute_Htot_basis_(Jc, Jb, mvp, i_freq);
                otherwise
                    error('Unknown coupling type! \n');
            end
           
            obj.task_solution_.set_H_total(obj.Htot);
            
        end
        
        % ------------------------------------------------------------- %
        
        function [] = compute_P_abs_(obj)
            % computes and stores the absorbed power
            
            dims.ql = obj.dims.ql;

            Pabs = obj.Etot .* conj(obj.task_solution_.Jb_plrz_tensor_);
            
            if obj.dims.ql == 12
                
                G = repmat([1; 1/12; 1/12; 1/12;], obj.dims.q,1);
                
                for i = dims.ql
                    Pabs(:,:,:,i) = G(i) .* Pabs(:,:,:,i);
                end
                
            end
            
            Pabs = sum(Pabs,4);
            Pabs = 0.5 * real(sum(Pabs(:))) * obj.scatterer.dom_vie.res.^3;
            
            obj.task_solution_.set_absorbed_power(Pabs)
        end
        
        % ------------------------------------------------------------- %
        function [SAR] = compute_SAR_(obj)
            
            rho     = obj.scatterer.prop_vie.rho;
            sigma_e = obj.scatterer.prop_vie.sigma_e;
            denom   = 2 * rho .* sigma_e; 
            Volume  = obj.scatterer.dom_vie.res.^3;
            
            idxsar = find(abs(denom(:))>1e-10);            
            Jb_tnsr = obj.task_solution_.Jb_plrz_tensor_;
            SAR.loc = Jb_tnsr .* conj(Jb_tnsr);
            
            % scale linear components (Gramian)
            if obj.dims.ql == 12
                G = repmat([1; 1/12; 1/12; 1/12;], obj.dims.q,1);
                for i = dims.ql
                    SAR.loc(:,:,:,i) = G(i) .* SAR.loc(:,:,:,i);
                end
            end
             
            SAR.loc = sum(SAR.loc,4);
            SAR.loc(idxsar) = SAR.loc(idxsar) ./ denom(idxsar);
            SAR.glob = sum(SAR.loc(:)) * Volume; 
            
            obj.task_solution_.set_SAR(SAR);
            
        end
        
        
        % ------------------------------------------------------------- %
        
        function [] = compute_b1p_(obj)
            % computes and stores b1+
            
            dims.ql = obj.dims.ql;
            
            mu0 = src_utils.EM_utils.mu0;
            
            B1_plus = mu0*abs(obj.Htot(:,:,:,1,:,:)+1i*obj.Htot(:,:,:,dims.ql/3+1,:,:));
            
            obj.task_solution_.set_B1_plus(B1_plus)
        end
        
        % ------------------------------------------------------------- %
        
        function [] = compute_b1m_(obj)
            % computes and stores b1-
            
            dims.ql = obj.dims.ql;
            
            mu0 = src_utils.EM_utils.mu0;
            
            B1_minus = mu0*abs(obj.Htot(:,:,:,1,:,:)-1i*obj.Htot(:,:,:,dims.ql/3+1,:,:));
            
            obj.task_solution_.set_B1_minus(B1_minus)
        end
    
    end
end