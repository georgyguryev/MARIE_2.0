classdef Assembly_VIE < src_assembly.Assembly_Base
% Class for body-specific problem assembly
% Author: Georgy Guryev, Cambridge, MA, 2019

    properties
        
        projector
        
        precorrection
        
        precond
                
        scatterer
        
        mvp

        
    end
    
% ====================================================================== %
    
methods
   
    function [obj] = Assembly_VIE(task_settings, scatterer, dims, coil)
        
        % call superclass constructor
        obj = obj@src_assembly.Assembly_Base(task_settings, dims);
        
        % get reference to body model
        obj.scatterer = scatterer;
        
        % get reference to coil model
        if ~isempty(coil)
            obj.coil = coil;
        end
        
        % instantiate preconditioner
        obj.precond = src_solver.Precond_LU(dims);
                
        % update sizer info
        obj.get_body_dimensions_();
        
        % get previous domain 
        prev_dom_vie = obj.setget_PrevDomain_vie();
        
        % set prev domain to current one if 
        if isempty(prev_dom_vie)
            obj.setget_PrevDomain_vie(obj.scatterer.dom_vie.r);
        end
        
    end
    
    % ------------------------------------------------------------- %

    function assemble_system(obj, freq)
        
        % compute material properties
        obj.compute_material_properties_vie_(freq);
        
        % assemble VIE and VSIE operators
        obj.assemble_operators_(freq);
        
        % assemble matrix-vector products
        obj.assemble_core_mvps_();
        
        % assemble precorrection matrices for implicit coupling
        obj.assemble_precorrection_vie_(freq);
        
        % assemble preconditioner for coupled/decompled solver
        obj.assemble_precond_(freq);
        
        % assemble vie-specific matrix-vector products
        obj.assemble_vie_mvps_(freq);
        
    end 
    
    % ------------------------------------------------------------- %
   
    function vie_solver = setup_solver(obj)
    % public method setup_solver() 

        vie_solver = src_solver.Solver_VIE(obj.mvp, obj.precond, obj.task_settings_);
            
        % setup relevant solver parameters
        vie_solver.setup_solver();
        
    end
    
    % ------------------------------------------------------------- %
    
    function [rhs_Nop, rhs_Kop] = generate_rhs(obj,freq)
        
        rhs_path     = fullfile('.','data', 'rhs');
        rhs_filename = obj.task_settings_.vie.Excitation_rhs;

        % function generates and saves rhs 
        [rhs_Nop, rhs_Kop] = obj.form_rhs(freq);
        
        % save rhs to mat file  
        save(fullfile(rhs_path,rhs_filename),'rhs_Nop', 'rhs_Kop');
    end
    
    % ------------------------------------------------------------- %
    
    function set_rhs(obj,freq)
        % function sets the rhs for VIE solver
        
        excitation_type = obj.task_settings_.vie.Excitation_type;
        
        rhs_path     = fullfile('.','data', 'rhs');
        rhs_filename = obj.task_settings_.vie.Excitation_rhs;
        rhs_fullname = fullfile(rhs_path, rhs_filename);
        
        switch excitation_type
        % If rhs exists in data/rhs -> load rhs; 
        % generate from scratch otherwise 
            case 'Plane_wave'
                
                % plarization and direction of propagation of a plane wave
                Eo = [1 0 0];
                k  = [0 0 1]; 
                
                [rhs_Nop, rhs_Kop] = obj.generate_plane_wave_excitation(Eo, k, freq);
                obj.dims.N_feeds = 1;
                
            case 'Coil'
                
                if exist(rhs_fullname, 'file')
                    load(rhs_fullname, 'rhs_Nop', 'rhs_Kop');
                else
                    [rhs_Nop, rhs_Kop] = obj.generate_rhs(freq);
                end
        end
        
        obj.rhs      = rhs_Nop;
        obj.rhs_Hinc = rhs_Kop;
    end
    

    
end
    
% ====================================================================== %
    
methods (Access = protected)
    
   % ------------------------------------------------------------- %
   function assemble_precorrection_(obj, freq)
        % method assemble_precorrection_() assembles the precorrection
        % matrices for implicit couplign; skips code otherwise
        
        % get coupling type
        coupling_type = obj.task_settings_.vsie.Coupling_Mode;
        
        % if implicit coupling is selected, compute precorrection mats.
        if strcmp(coupling_type, 'pFFT')
            
            % instantiate precorrection matrices
            obj.precorrection = src_precorrection.Precorrection(obj.coil, obj.scatterer,...
                obj.projector, obj.task_settings_,...
                obj.dims, obj.mvp, obj.Zcoil_, freq);
            
            % compute Coil2Body precorrection
            obj.precorrection.assemble_precorrection();
            
        end
    end

    
    % ------------------------------------------------------------- %

    function solver = setup_iterative_solver_(obj, mvp)
        % private method setup_solver()
        
        % specify solver type
        solver_type = obj.task_settings_.vie.Iterative_method;
        
        % setup appropriate solver
        switch solver_type
            case 'GMRES_DR'
                solver = src_solver.gmres_DR_solver(@(Jin) mvp(Jin), obj.precond, obj.task_settings_);
            case 'GMRES'
                solver = src_solver.gmres_solver(@(Jin) mvp(Jin), obj.precond, obj.task_settings_);
            case 'BiCGSTAB'
                solver = src_solver.bicgstab_solver(@(Jin) mvp(Jin), obj.precond, obj.task_settings_);
        end
        
    end
    
    
    % ------------------------------------------------------------- %
    
    function assemble_operators_(obj, freq)
        
        excitation_type = obj.task_settings_.vie.Excitation_type;
        rhs_filename    = obj.task_settings_.vie.Excitation_rhs;
        coupling_mode   = obj.task_settings_.vsie.Coupling_Mode;
        
        assemble_rhs_flag = (strcmp(excitation_type, 'Coil')) && ...
                            (~exist(rhs_filename));
        
        % initialize projector matrix
        obj.projector = src_projector.Projector();
        
        % compute material properties for vie
        obj.compute_material_properties_vie_(freq);
        
        % assemble operators for implicit nested VSIE
        obj.assemble_vie_operators_(freq);
        
        % get operator dimensions
        obj.get_operator_dimensions_vie_();
        
        % assemble coupling operators if there is no rhs
        if assemble_rhs_flag
            
            % call some assebly method
            [obj.Zcoil_, obj.Fcoil_, obj.feed_tune_ports_] = Assembly_SIE_par(obj.coil, freq);

            switch coupling_mode
                case 'pFFT'
                    
                    obj.assemble_operators_vie_coupling_pFFT_(freq);
                    
                otherwise

                    obj.assemble_operators_vie_coupling_Dense_(freq);
            end
            
        end
        
    end
    
    % ------------------------------------------------------------- %
    
    function assemble_vie_operators_(obj, freq)
        % method assemble_implicit_nested_() assembles operators
        % and coupling matrices for nested implicit coupling
        
        % get gpu flag
        gpu_flag    = obj.task_settings_.general.GPU_flag;
        tucker_flag = obj.task_settings_.vie.Tucker;
        
        obj.operators.N_vie_ = obj.setget_Prev_N_vie();

        % generate N operator for VIE domain
        if ~isequal(obj.setget_PrevDomain_vie(), obj.scatterer.dom_vie.r) || ...
                    isempty(obj.operators.N_vie_) ||...
                    (obj.setgetPrevFreq_VIE() ~= freq) || ...
                    isempty(obj.operators.K_vie_)
            
            % assemble operators for the VIE domain
            obj.operators.N_vie_ = obj.assemble_Nop_(obj.scatterer.dom_vie, freq);
            obj.operators.K_vie_ = obj.assemble_Kop_(obj.scatterer.dom_vie, freq);
            
            % save the operator to persistent variable
            obj.setget_Prev_N_vie(obj.operators.N_vie_);
            obj.setget_PrevDomain_vie(obj.scatterer.dom_vie.r);
        end

        obj.setgetPrevFreq_VIE(freq);

        % get operator dimensions
        obj.get_operator_dimensions_vsie_();
        
        % if gpu flag is nonzero allocate gpu device
        if gpu_flag
            
            % get a handle to GPU
            gpuDevice(gpu_flag);
            
            % allocate operators N operators for vie and extended domain on GPU
            obj.operators.N_vie_gpu  = src_utils.to_GPU(obj.operators.N_vie_, tucker_flag);
            
            % move material properties to GPU
            obj.scatterer.prop_vie.Mcr = src_utils.to_GPU(obj.scatterer.prop_vie.Mcr,0);
            obj.scatterer.prop_vie.Mrc = src_utils.to_GPU(obj.scatterer.prop_vie.Mrc,0);
            
        end
        
    end
    
    % ------------------------------------------------------------- %
    
    function assemble_core_mvps_(obj)
        % method assemble_mvps_() calls vie-specific  constructor to set
        % the appropriate mvps for vie task
        
        % construct mvp_vie() 
        obj.mvp = src_mvp.mvp_vie(obj.task_settings_, obj.operators,...
                                  obj.scatterer, obj.projector,...
                                  obj.Zcoil_,obj.dims);
                              
    end
    
    % ------------------------------------------------------------- %
    
    function assemble_vie_mvps_(obj, freq)
        % method assemble_vie_mvps_() sets up vie mvps for VIE solver
        
        % set up mvp for vie solver
        obj.mvp.set_final_mvp_vie(obj.precorrection, freq);
    end
    
    % ------------------------------------------------------------- %
    
    function [Mr, Mc, Mr_inv, Mc_inv, Mcr_inv] = compute_properties_matrix_(obj, mat_prop, index,freq)
        
        % compute omega
        emu = src_utils.EM_utils(freq);
                
        % define complex relative permittivity and susceptibility
        Mr  = mat_prop.epsilon_r + mat_prop.sigma_e / (emu.ce);
        Mc  = (Mr - 1);
        Mcr = Mc ./ Mr; 
        
        % allocate size for matrix of inverse susceptibilities
        Mr_inv  = zeros(size(Mr));
        Mc_inv  = zeros(size(Mc));
        Mcr_inv = zeros(size(Mcr));
        
        % compute inverse of the suseptibility matrix
        Mr_inv(index.S_1d)  = 1 ./ Mr(index.S_1d);
        Mc_inv(index.S_1d)  = 1 ./ Mc(index.S_1d);
        Mcr_inv(index.S_1d) = 1 ./ Mcr(index.S_1d); 
        
    end
    
    % ------------------------------------------------------------- %
    
    function compute_material_properties_vie_(obj, freq)
        
        % compute matrices of material properties for vie domain
        [Mr_vie, Mc_vie, Mr_inv_vie, Mc_inv_vie, Mcr_inv_vie] = obj.compute_properties_matrix_(...
                                                       obj.scatterer.prop_vie,...
                                                       obj.scatterer.index_vie,freq);
        % store material properties                                            
        obj.scatterer.prop_vie.set_material_matrices(Mr_vie, Mc_vie, Mr_inv_vie, Mc_inv_vie, Mcr_inv_vie);                                           
    end
    
    % ------------------------------------------------------------- %
    
    function N_op = assemble_Nop_(obj, dom, freq)
        % method assemble_Nop_(freq) assembles N operator
        
        tucker_flag = obj.task_settings_.vie.Tucker;
        tucker_tol  = obj.task_settings_.vie.TuckerTolerance;

        N_op = src_operator.getOPERATORS(dom, freq, 'N', obj.dims.l, tucker_flag, tucker_tol);
   
    end
    
    % ------------------------------------------------------------- %
    
    function K_op = assemble_Kop_(obj, dom, freq)
        % method assemble_Kop_(freq) assembles K operator
        
        tucker_flag = obj.task_settings_.vie.Tucker;
        tucker_tol  = obj.task_settings_.vie.TuckerTolerance;
        
        K_op = src_operator.getOPERATORS(dom, freq, 'K', obj.dims.l, tucker_flag, tucker_tol);
    end
    
    % ------------------------------------------------------------- %
    
    function get_body_dimensions_(obj)
        % private function updates obj.dims fields that correspond to body
        % dimensions
        obj.dims.get_body_dims(obj.scatterer.dom_vie.r);
            
    end
    
    % ------------------------------------------------------------- %
    
    function [op_dims] = get_operator_dimensions_(obj, operator)
        % update operator-related dimensions
        
         % update obj.dims components        
        if obj.task_settings_.vie.Tucker
            
            L = size(operator(1,1).U1,1);
            M = size(operator(1,1).U2,1);
            N = size(operator(1,1).U3,1);

        else
            [L,M,N,~] = size(operator);
        end
        
        % return dimensions
        op_dims = [L,M,N];
    end
    
    % ------------------------------------------------------------- %
    
    function get_operator_dimensions_vie_(obj)
        % get operator dimensions
        
        if ~isempty(obj.operators.N_vie_)
            obj.dims.op_vie = obj.get_operator_dimensions_(obj.operators.N_vie_);
        else
            obj.dims.op_vie = obj.get_operator_dimensions_(obj.operators.K_vie_);
        end
        
    end
    
    % ------------------------------------------------------------- %
    
    function get_near_dimensions_(obj)
        
        
        L_max = max(max(obj.coil.Ln));
        N_exp = ceil(2 * L_max / obj.scatterer.dom_vie.res); 
        
        error_str = 'The coil resolution is too coarse for pFFT, please refine mesh or choose another method';
        
         if obj.dims.ql == 3
            if N_exp <= 3
                N_exp = 3;
            elseif N_exp <= 5
                N_exp = 5;
            else 
                error(error_str);
            end
        else 
            if N_exp <= 3
                N_exp = 3;
            else 
                error(error_str);
            end
        end
        
        % diameter of a bounding box
        N_near = 6 * ceil(N_exp/2) + 1;
        
        % get 1D dimensions of expansion and near domains
        obj.dims.N_exp_1D  = N_exp;
        obj.dims.N_near_1D = N_near;
        
        % get 3D dimensions of expansion and near domains
        obj.dims.N_exp_3D  = N_exp.^3;
        obj.dims.N_near_3D = N_near.^3; 
        
%         obj.projector.near_boundary_width = obj.task_settings_.vsie.Near_boundary_width;
        
        %
        obj.dims.near = [obj.dims.N_near_1D, obj.dims.N_near_1D, obj.dims.N_near_1D, obj.dims.ql];
        obj.dims.exp  = [obj.dims.N_exp_1D, obj.dims.N_exp_1D, obj.dims.N_exp_1D, obj.dims.ql];
    end
    
    % ------------------------------------------------------------- %
    
    function extend_domain_(obj, freq)
                
        % get dimensions for near domain
        obj.get_near_dimensions_()
        
        % function extends vie to ext domain and resets all properties
        [dom_ext, prop_ext] = src_scatterer.extend_domain(obj.coil, obj.scatterer,...
            obj.projector, obj.dims,...
            freq);
        
        % set extended domain
        obj.scatterer.set_extended_domain(dom_ext, prop_ext);
        
        % get previous domain
        if isempty(obj.setget_PrevDomain_ext())
            obj.setget_PrevDomain_ext(dom_ext.r);
        end
        
        
        % form near lists
        obj.projector.form_near_lists(obj.coil, obj.scatterer, obj.projector, obj.dims);
               
        % generate extended domain
        obj.scatterer.generate_near_domain();
        
    end
    
    % ------------------------------------------------------------- %
    
    function assemble_projection_matrix_(obj, freq)
        % function projector_assembly() serves as interface between various modes
        % of projector assembly. Function returns assembled projection matrix PS
        % for the specified collocation mode
        
        obj.projector.assemble_projection_matrix(obj, freq);
        
    end
    
    % ---------------------------------------------------------------------- %
    
    function assemble_coupling_matrix_(obj, freq)
        
        %
        gpu_flag = obj.task_settings_.basis.GPU_flag;
        
        [Zbc_Nop, Zbc_Kop] = src_coupling.assemble_coupling_matrices(obj.task_settings_, obj.scatterer, obj.coil, ...
            obj.dims, freq, [], [], 0);
        if gpu_flag
            
            obj.operators.Zbc_Nop = gpuArray(Zbc_Nop);
            obj.operators.Zbc_Kop = gpuArray(Zbc_Kop);
        else
            obj.operators.Zbc_Nop = Zbc_Nop;
            obj.operators.Zbc_Kop = Zbc_Kop;
        end
        
        clear Zbc_Nop Zbc_Kop;
    end
    
end

% ====================================================================== %

methods(Access = private)
    
    % ------------------------------------------------------------------- %
    
    function get_operator_dimensions_vsie_(obj)
        
        % Get vsie-related operator dimensions
        obj.get_operator_dimensions_vie_();

    end
    
    % ------------------------------------------------------------- %
    function assemble_precorrection_vie_(obj, freq)
        % method assemble_precorrection_() assembles the precorrection
        % matrices for implicit couplign; skips code otherwise
        
        % get coupling type
        excitation_type = obj.task_settings_.vie.Excitation_type;
        rhs_filename    = obj.task_settings_.vie.Excitation_rhs;
        
        assemble_rhs_flag = (strcmp(excitation_type, 'Coil')) && ...
                            (~exist(rhs_filename));
        
        % if implicit coupling is selected, compute precorrection mats.
        if assemble_rhs_flag
            
            obj.assemble_precorrection_(freq);
        end
    end
    
    % ------------------------------------------------------------------- %
    
    function assemble_precond_(obj, freq)
        % mathod assemble_precond_() assembles preconditioner for
        % coupled or decoupled systems
        
        idx_1d      = obj.scatterer.index_vie.S_1d;
        Mcr_vie_inv = obj.scatterer.prop_vie.Mcr_inv;
        
        % VIE preconditioner
        obj.precond.assemble_precond_vie(idx_1d, Mcr_vie_inv, obj.mvp, freq);
                
    end
    
    % ------------------------------------------------------------------- %
    
    function [rhs_Nop, rhs_Kop] = form_rhs(obj, freq)
        % method form_rhs()
        
        
        N_vie = obj.dims.N_scat  * obj.dims.ql;
        rhs_Nop = zeros(N_vie, obj.dims.N_feeds);
        rhs_Kop = zeros(N_vie, obj.dims.N_feeds);

        % call some assebly method
        [obj.Zcoil_, obj.Fcoil_, obj.feed_tune_ports_] = Assembly_SIE_par(obj.coil, freq);
                
        % compute coil currents without a scatterer
        Jc_free_space = obj.Zcoil_ \ obj.Fcoil_;
        
        for port = 1:obj.dims.N_feeds
            % compute Einc due to surface currents
            rhs_Nop_port = obj.mvp.c2b_coupling_vie_Nop(Jc_free_space(:,port));
            rhs_Kop_port = obj.mvp.c2b_coupling_vie_Kop(Jc_free_space(:,port));
            
            rhs_Nop(:,port) = rhs_Nop_port;
            rhs_Kop(:,port) = rhs_Kop_port;
        end

    end
    
    
    % ------------------------------------------------------------------- %
    function [E_inc, H_inc] = generate_plane_wave_excitation(obj, Eo, k, freq)
        
        % get em constants
        emu = src_utils.EM_utils(freq);
        
        % get index 
        idxS_3D = obj.scatterer.index_vie.S_3d;
        
        N_vox = obj.dims.N_scat_vox;
        
        % define polarization and direction of propagation       
        Eo = Eo / norm(Eo);        
        k  = k / norm(k) * emu.k0;
        
        % generate fields         
        [E, H] = src_utils.PlaneWave_Excitation(obj.scatterer, k, emu.omega,...
                                                obj.dims, Eo);
                                                   
        if obj.dims.ql == 3                                   
            % scale by volume and EM scaling factor
            E_inc = E(idxS_3D); 
            H_inc = H(idxS_3D);
        elseif obj.dims.ql == 12
            rhs_E = E(idxS_3D);
            rhs_H = H(idxS_3D);
            
            E_inc = zeros(obj.dims.l * size(rhs_E,1),1);
            H_inc = zeros(obj.dims.l * size(rhs_H,1),1);  
            
            E_inc(1:N_vox) = rhs_E(1:N_vox);
            H_inc(1:N_vox) = rhs_H(1:N_vox);
            E_inc(4*N_vox+1:5*N_vox)  = rhs_E(N_vox + 1:2*N_vox);
            H_inc(4*N_vox+1:5*N_vox)  = rhs_H(N_vox + 1:2*N_vox);
            E_inc(8*N_vox+1:9*N_vox)  = rhs_E(2*N_vox + 1:3*N_vox);
            H_inc(8*N_vox+1:9*N_vox)  = rhs_H(2*N_vox + 1:3*N_vox);
         end
        
    end
    
     % ------------------------------------------------------------------- %
   
    function [] = assemble_operators_vie_coupling_pFFT_(obj,freq)
        
        gpu_flag    = obj.task_settings_.general.GPU_flag;
        tucker_flag = obj.task_settings_.vie.Tucker;
        
        obj.extend_domain_(freq);        
        
        % assemble operators for the Ext domain
        obj.operators.N_ext_ = obj.assemble_Nop_(obj.scatterer.dom_ext, freq);
        obj.operators.K_ext_ = obj.assemble_Kop_(obj.scatterer.dom_ext, freq);
        
        % assemble operators for the Near domain
        obj.operators.N_near_ = obj.assemble_Nop_(obj.scatterer.dom_near, freq);
        obj.operators.K_near_ = obj.assemble_Kop_(obj.scatterer.dom_near, freq);
        
        % assemble projection matrices
        obj.assemble_projection_matrix_(freq);
        
        % get operator dimensions
        obj.dims.op_ext  = obj.get_operator_dimensions_(obj.operators.N_ext_);
        obj.dims.op_near = obj.get_operator_dimensions_(obj.operators.N_near_);
        
        % if gpu flag is nonzero allocate gpu device
        if gpu_flag
            
            % allocate operators N operators for vie and extended domain on GPU
            obj.operators.N_ext_gpu  = src_utils.to_GPU(obj.operators.N_ext_, tucker_flag);
            obj.operators.K_ext_gpu  = src_utils.to_GPU(obj.operators.K_ext_, tucker_flag);
            
            % obj.operators.N_vie_gpu  = src_utils.to_GPU(obj.operators.N_vie_, tucker_flag);
            obj.operators.N_near_gpu = src_utils.to_GPU(obj.operators.N_near_, tucker_flag);
            obj.operators.K_near_gpu = src_utils.to_GPU(obj.operators.K_near_, tucker_flag);
            
        end
        
    end

    % ------------------------------------------------------------------- %

    function [] = assemble_operators_vie_coupling_Dense_(obj, freq)
        obj.assemble_coupling_matrix_(freq);
    end
               
end
 
   
end