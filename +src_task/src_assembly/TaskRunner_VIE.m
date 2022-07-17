classdef TaskRunner_VIE < TaskRunner_Base
    
    properties
                
        % store body model
        coil
        scatterer
        
        % assembly_vie object
        assembly_vie        
        
        % iterative solver
        solver 
        
        % temporary var. for resulting currents
        Jcb_sol
        
        % 
        sol_time
        rel_res
        res_vec = {}
        
        
    end
    
    methods
        
        % -------------------------------------------------------------- %
        
        % define Task_VIE constructor
        function [obj] = TaskRunner_VIE(task_settings)
            
            % call base class constructor
            obj = obj@TaskRunner_Base(task_settings);
            
            if strcmp(obj.task_settings_.vie.Excitation_type, 'Coil')
            
                % set coil loader and load a coil (for Excitation)
                obj.coil_loader = src_geo.Coil_Loader(task_settings.sie);
                
                obj.load_coil_model_();
            end

            % load coil and body models            
            obj.load_scatterer_(); 
            
            obj.assembly_vie = src_assembly.Assembly_VIE(task_settings, obj.scatterer,...
                                                         obj.dims, obj.coil);

        end

        % -------------------------------------------------------------- %

        % solve the problem
        function Jsol = run(obj)
            
            % get working frequencies
            [freqs, N_freq] = obj.get_working_frequencies_();
     
            % run frequency sweep 
            for  n_freq = 1:N_freq
                
                % get current frequency
                freq = freqs(n_freq);
                
                %assemble VIE system
                obj.assemble_VIE_system_(freq);
                
                % setup iterative solver
                obj.solver = obj.assembly_vie.setup_solver();
                
                if isempty(obj.assembly_vie.rhs)
                    % generate rhs
                   obj.assembly_vie.set_rhs(freq);
                end
                
                for feed_port = 1:obj.dims.N_feeds
                    
                    fprintf('Solving VIE problem for Port %d \n', feed_port);

                    tic
                    
                    current_rhs = obj.assembly_vie.rhs(:,feed_port);
                    [Jcb, relative_res, res_vector] = obj.solver.run(current_rhs, []);
                    
                    solve_time = toc;
                    
                    % store results in temporary tensors
                    obj.Jcb_sol(:,feed_port, n_freq) = Jcb;
                    obj.sol_time(feed_port, n_freq)  = solve_time;
                    obj.rel_res(:,feed_port, n_freq) = relative_res;
                    obj.res_vec(feed_port, n_freq)   = {res_vector};
                    
                    fprintf('VIE problem solved in %d iterations\n', length(res_vector));
                    fprintf('Relative residual: %e \n', relative_res);
                    fprintf('Elapsed time: %4.2f seconds \n\n', solve_time);

                end
                                                
                % compute electric, magnetic fields
                fprintf('Computing Electric field ... \n\n');

                obj.compute_E_total_(n_freq);
                
                fprintf('Computing Magnetic field ... \n\n');
                
                obj.compute_H_total_VIE_(n_freq);
            end
            
            % update solution object
            obj.update_solution_();
            
            obj.compute_P_abs_();
            obj.compute_SAR_();

        end

    end
    
    
    methods (Access = protected)
        
        % -------------------------------------------------------------- %
       
        function assemble_VIE_system_(obj, freq)
            obj.assembly_vie.assemble_system(freq);
        end
        
        % -------------------------------------------------------------- %
        
        function load_scatterer_(obj)
            % function load_scatterer(obj, freq) 
            
            % get current frequency
            [freqs, ~] = obj.get_working_frequencies_();
                        
            % load RHBM structure from file
            RHBM = load(obj.task_settings_.vie.BodyModel);
            
            % init scatterer object
            obj.scatterer = src_scatterer.Scatterer(RHBM, obj.dims, freqs(1));
            
        end
        
        % -------------------------------------------------------------- %
        
    end
    
    methods (Access = private)
        
        % -------------------------------------------------------------- %
        
        function [] = load_coil_model_(obj)
            % private method loads coil to task runner
            obj.coil = obj.coil_loader.load_coil();
            
            % get coil dimensions
            obj.dims.get_coil_dims(obj.coil);
        end
        
        % -------------------------------------------------------------- %
        
        function update_solution_(obj)
            
            % get dimensions of the problem
            dims = obj.dims;
            
            [freqs, N_freq] = obj.get_working_frequencies_();
            idx_ql   = obj.scatterer.index_vie.index_ql(dims.ql, dims.Nvox_vie);
            
            Jb_vec    = obj.Jcb_sol(:,:,:);
            Jb_tensor = zeros([dims.vie, dims.N_feeds, N_freq]);

            for i = 1:1%dims.N_feeds
                for j = 1:N_freq
                    
                    N_shift = dims.Nvox_vie * dims.ql;
                    idx = idx_ql + (i-1) * N_shift + ...
                          (j-1) * N_shift * dims.N_feeds;
                    Jb_tensor(idx) = Jb_vec(:,i,j);
                end
            end    
            
            % store solution of the VSIE problem in task_solution
            obj.task_solution_.set_plrz_currents_vec(Jb_vec);
            obj.task_solution_.set_plrz_currents_tensor(Jb_tensor);
            
            % store relative residual and res. vector
            obj.task_solution_.set_relative_residual(obj.rel_res);
            obj.task_solution_.set_residual_vector(obj.res_vec);
            obj.task_solution_.set_lapse_time(obj.sol_time);
            
            % store total electric and magnetic fields
            obj.task_solution_.set_E_total(obj.Etot)
            obj.task_solution_.set_H_total(obj.Htot)

            % set the relative permittivity and conductivity
            obj.task_solution_.set_epsilon_r(obj.scatterer.prop_vie.epsilon_r);
            obj.task_solution_.set_sigma_e(obj.scatterer.prop_vie.sigma_e);
            
            % set the frequencies
            obj.task_solution_.set_freqs(freqs);
            
        end
        
        % -------------------------------------------------------------- %
        
        function [] = compute_E_total_(obj, i_freq)
        % compute total electric fields 
        
            % get tissue currents 
            obj.compute_E_total_base_(obj.Jcb_sol, obj.assembly_vie.mvp,...
                                       i_freq)
        
        end

        % -------------------------------------------------------------- %
        
        function [] = compute_H_total_VIE_(obj, i_freq)
            % computes the magnetic field
            
            % get dimensions of the problem
            dims = obj.dims;
            n1 = dims.L_vie;
            n2 = dims.M_vie;
            n3 = dims.N_vie;
            n_feeds = dims.N_feeds;
            
            res = obj.scatterer.dom_vie.res;
            
            % compute tested total electric field
            t_Htotal = zeros(n1,n2,n3,dims.ql);
                        
            idxS = obj.assembly_vie.scatterer.index_vie.index_ql(dims.ql, dims.Nvox_vie);
            gpu_flag = obj.task_settings_.general.GPU_flag;
                                   
            for feed_port = 1:dims.N_feeds
                                
                % get the polarization currents
                Jsol         = zeros(n1 * n2 * n3 * dims.ql, n_feeds);
                Jsol(idxS,:) = obj.Jcb_sol(:,:, i_freq);
                Jsol         = reshape(Jsol,[n1, n2, n3, dims.ql, n_feeds]);
                
                % compute tested incident magnetic field
                t_Hinc = obj.assembly_vie.rhs_Hinc(:,feed_port);
     
                % compute tested scattered magnetic field
                t_Hsca = obj.assembly_vie.mvp.K_vie(Jsol(:,:,:,:,feed_port));
                t_Hsca = t_Hsca(idxS);
                
                % compute total tested magnetic field
                if gpu_flag
                    t_Htotal(idxS) = gather(t_Hsca + t_Hinc);
                else
                    t_Htotal(idxS) = t_Hsca + t_Hinc;
                end
                
                % PWL scaling
                if dims.ql == 12
                    
                    G = repmat([1; 1/12; 1/12; 1/12;], obj.dims.q,1);
                    
                    for i = dims.ql
                        t_Htotal(:,:,:,i) = G(i).*t_Htotal(:,:,:,i);
                    end
                end
                
                % store only the total magnetic field
                obj.Htot(:,:,:,:,feed_port,i_freq) = t_Htotal / (res^3);
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
                for i = obj.dims.ql
                    SAR.loc(:,:,:,i) = G(i) .* SAR.loc(:,:,:,i);
                end
            end
             
            SAR.loc = sum(SAR.loc,4);
            SAR.loc(idxsar) = SAR.loc(idxsar) ./ denom(idxsar);
            SAR.glob = sum(SAR.loc(:)) * Volume; 
            
            obj.task_solution_.set_SAR(SAR);
            
        end
        
    end    
    
end