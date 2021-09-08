classdef Assembly_VSIE < src_assembly.Assembly_VIE & src_assembly.Assembly_SIE
% Class for full EM problem assembly (coil + body)
% Author: Georgy Guryev, Cambridge, MA, 2019    
    properties
        
        coupling

    end
    
% ====================================================================== %
    
    methods
        
        function obj = Assembly_VSIE(task_settings, coil, scatterer, dims)
            
            % call Assembly SIE and Assembly VIE contructors
            obj = obj@src_assembly.Assembly_SIE(task_settings, coil, dims);
            obj = obj@src_assembly.Assembly_VIE(task_settings, scatterer, dims, coil);
            
            % instantiate preconditioner 
            obj.precond = src_solver.Precond_LU(dims);
            
        end
        
        % --------------------------------------------------------------- %

        function assemble_system(obj, freq)

            % assemble edge2edge SIE matrix 
            obj.assemble_SIE_system_(freq);
                        
            % assemble VIE and VSIE operators
            obj.assemble_operators_(freq);  
            
            % assemble matrix-vector products
            obj.assemble_core_mvps_();
                      
            % assemble precorrection matrices for implicit coupling
            obj.assemble_precorrection_(freq);
            
            % assemble preconditioner for coupled/decompled solver
            obj.assemble_precond_(freq);
            
            % assemble final mvp for solver
            obj.assemble_final_mvp_(freq);
            
            % form rhs for iterative solver
            obj.form_rhs(freq);
            
        end
        
        % --------------------------------------------------------------- %

        function vsie_solver = setup_solver(obj)
        % public method setup_solver() 
                           
            % get formulation type
            formulation_type = obj.task_settings_.vsie.Solver_mode;
            
            % construct appropriate solver type 
            switch formulation_type
                
                case 'Coupled'
                    vsie_solver = src_solver.Solver_Coupled(obj.mvp, obj.precond, obj.task_settings_);
                case 'Decoupled'
                    vsie_solver = src_solver.Solver_Decoupled(obj.mvp, obj.precond, obj.task_settings_);
            end
            
            % setup relevant solver parameters
            vsie_solver.setup_solver();
        end
        
        % --------------------------------------------------------------- %

        function fP_mvp_near = assemble_near_mvp(obj)
            
            % get N operator for near domain
            N_near = obj.operators.N_near_;
            
            % set size of N operator for near domain
            if obj.task_settings_.vie.Tucker
                n1 = size(N_near(1,1).U1,1);
                n2 = size(N_near(1,1).U2,1);
                n3 = size(N_near(1,1).U3,1); 
                obj.dims.Nop_near = [n1 n2 n3];
            else
                obj.dims.Nop_near = size(N_near);
            end
            
            % get
            if obj.task_settings_.vie.PWX
                fN_mvp_near = @(Jin) src_mvp.N_mvp_pwl(Jin, N_near, obj.dims.near, obj.dims.Nop_near, obj.dims);
            else
                fN_mvp_near = @(Jin) src_mvp.N_mvp_pwc(Jin, N_near, obj.dims.near, obj.dims.Nop_near);
            end
            
            fP_mvp_near = @(Jin) src_mvp.mvp_P(Jin, fN_mvp_near, obj.projector, obj.dims);
        end
        
    end
    
 % ====================================================================== %
   
    methods (Access = protected)
        
        function assemble_operators_(obj, freq)
                        
            % get solver and coupling modes
            coupling_type  = obj.task_settings_.vsie.Coupling_Mode;
            
            % add contribution of external circuits to Zcc
            obj.add_matching_impedance_(freq);
            
            % form rhs taking into account ext. circuitry
            obj.form_rhs_sie_();
                        
            obj.Zcoil_inv_ = eye(size(obj.Zcoil_)) / obj.Zcoil_; 
            
            obj.operators.Jc_ini =   obj.Zcoil_inv_ * obj.Fcoil_;

            if strcmp(coupling_type ,'Dense')
                
                % compute material properties for vie
                obj.compute_material_properties_vie_(freq);
                
                % assemble N and K operators for explicit coupling
                obj.assemble_operators_explicit_(freq);                    
                
            elseif strcmp(coupling_type, 'pFFT')
                % assemble N and K operators for implicit decoupled
                
                tic; 
                % initialize projector matrix
                obj.projector = src_projector.Projector();
                toc;
                
                % extend domain
                obj.extend_domain_(freq);
                
                % assemble material properties for extended domain
                obj.compute_material_properties_ext_(freq);

                % compute material properties for vie
                obj.compute_material_properties_vie_(freq);
              
                % assemble operators for implicit nested VSIE
                obj.assemble_operators_pFFT_(freq);
                
            elseif strcmp(coupling_type, 'Basis')
                
                % compute material properties for vie
                obj.compute_material_properties_vie_(freq);
                
                obj.assemble_operators_basis_(freq);
            end
            
        end
        
        % --------------------------------------------------------------- %
        
        function assemble_core_mvps_(obj)
            % method assemble_mvps_() calls vie-specific  constructor to set
            % the appropriate mvps for vie task
            
            % construct mvp_vie()
            obj.mvp = src_mvp.mvp_vsie(obj.task_settings_, obj.operators,...
                                  obj.scatterer, obj.projector, obj.Zcoil_,...
                                  obj.Zcoil_inv_, obj.dims);
        end
        
    end
    
 % ====================================================================== %
    
    methods (Access = private)
        
        
        function assemble_precond_(obj, freq)
            % mathod assemble_precond_() assembles preconditioner for
            % coupled or decoupled systems

            % get solver mode
            coupling_method = obj.task_settings_.vsie.Coupling_Mode;
            solver_mode = obj.task_settings_.vsie.Solver_mode;
            
            if strcmp(solver_mode, 'Coupled')
                
                obj.precond.assemble_precond_sie(obj.Zcoil_);
                
                if strcmp(coupling_method, 'pFFT')
                    
                    idx_1d    = obj.scatterer.index_ext.S_1d;
                    Mcr_inv = obj.scatterer.prop_ext.Mcr_inv;
                else 
                    idx_1d    = obj.scatterer.index_vie.S_1d;
                    Mcr_inv = obj.scatterer.prop_vie.Mcr_inv;
                end
                
                obj.precond.assemble_precond_vsie_coupled(idx_1d, Mcr_inv, obj.mvp, freq);
                
            else
                
                idx_1d    = obj.scatterer.index_vie.S_1d;
                Mcr_vie_inv = obj.scatterer.prop_vie.Mcr_inv;
                
                % preconditioner for vsie (vie domain)
                obj.precond.assemble_precond_vie(idx_1d, Mcr_vie_inv, obj.mvp, freq);
            end
            
        end
        
        % ---------------------------------------------------------------------- %
        
        function compute_material_properties_ext_(obj, freq)
            % assemble material properties for extended domain
            [Mr_ext, Mc_ext, Mr_ext_inv, Mc_ext_inv, Mcr_ext_inv] = obj.compute_properties_matrix_(...
                                                           obj.scatterer.prop_ext,...
                                                           obj.scatterer.index_ext,freq);
                                                       
            % set matrices of material properties                                           
            obj.scatterer.prop_ext.set_material_matrices(Mr_ext, Mc_ext, Mr_ext_inv, Mc_ext_inv, Mcr_ext_inv);                                           
                                           
        end
        
        % ---------------------------------------------------------------------- %        
        
        function assemble_operators_explicit_(obj, freq)
            % method assemble_explicit_() assembles operators and coupling
            % matrices for nested/unnested explicit coupling
            
            % get gpu and tucker flags
            gpu_flag     = obj.task_settings_.general.GPU_flag;
            tucker_flag = obj.task_settings_.vie.Tucker;
            
            obj.operators.N_vie_  = obj.setget_Prev_N_vie();
            
            % generate N operator for VIE domain       
            if ~isequal(obj.setget_PrevDomain_vie(), obj.scatterer.dom_vie.r) || ...
                isempty(obj.operators.N_vie_) || (obj.setgetPrevFreq_VIE() ~= freq) || ...
                isempty(obj.operators.K_vie_)     
                
                % assemble operators for the VIE domain
                obj.operators.N_vie_ = obj.assemble_Nop_(obj.scatterer.dom_vie, freq);

                % save the operator to persistent variable
                obj.setget_Prev_N_vie(obj.operators.N_vie_);
                obj.setget_PrevDomain_vie(obj.scatterer.dom_vie.r);
                obj.setgetPrevFreq_VIE(freq);
                obj.operators.K_vie_ = obj.assemble_Kop_(obj.scatterer.dom_vie, freq);
            end
        
            % get operator dimensions
            obj.get_operator_dimensions_vie_();
            
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
            
            % assemble coupling matrix (discritization of coupling operator)
            obj.assemble_coupling_matrix_(freq);
            
        end
        
        % ---------------------------------------------------------------------- %
        
        function assemble_operators_basis_(obj, freq)
            % method assemble_implicit_nested_() assembles operators
            % and coupling matrices for nested implicit coupling
            
            % get gpu and tucker flags
            gpu_flag     = obj.task_settings_.general.GPU_flag;
            tucker_flag = obj.task_settings_.vie.Tucker;
            
            obj.operators.N_vie_  = obj.setget_Prev_N_vie();
            

            % generate N operator for VIE domain       
            if ~isequal(obj.setget_PrevDomain_vie(), obj.scatterer.dom_vie.r) || ...
                isempty(obj.operators.N_vie_) || (obj.setgetPrevFreq_VIE() ~= freq) || ...
                isempty(obj.operators.K_vie_)     
                
                % assemble operators for the VIE domain
                obj.operators.N_vie_ = obj.assemble_Nop_(obj.scatterer.dom_vie, freq);
    
                % save the operator to persistent variable
                obj.setget_Prev_N_vie(obj.operators.N_vie_);
                obj.setget_PrevDomain_vie(obj.scatterer.dom_vie.r);
                obj.setgetPrevFreq_VIE(freq);
                obj.operators.K_vie_ = obj.assemble_Kop_(obj.scatterer.dom_vie, freq);
                
            end
        
            % get operator dimensions
            obj.get_operator_dimensions_vie_();
            
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
            
            % 
            obj.assemble_basis_coupling_(freq);
        end
        
        % ---------------------------------------------------------------------- %
        
        function assemble_operators_pFFT_(obj, freq)           
            
            % get gpu flag
            gpu_flag    = obj.task_settings_.general.GPU_flag;
            tucker_flag = obj.task_settings_.vie.Tucker;

            % get operator from the previous run (if exists)
            obj.operators.N_vie_  = obj.setget_Prev_N_vie();
            obj.operators.N_ext_  = obj.setget_Prev_N_ext();
            
            
            % generate N operator for VIE domain       
            if ~isequal(obj.setget_PrevDomain_vie(), obj.scatterer.dom_vie.r) || ...
                isempty(obj.operators.N_vie_) || (obj.setgetPrevFreq_VIE() ~= freq) || ...
                isempty(obj.operators.K_vie_)
                
                
                % assemble operators for the VIE domain
                obj.operators.N_vie_ = obj.assemble_Nop_(obj.scatterer.dom_vie, freq);

                % save the operator to persistent variable
                obj.setget_Prev_N_vie(obj.operators.N_vie_);
                obj.setget_PrevDomain_vie(obj.scatterer.dom_vie.r);
                
                obj.operators.K_vie_ = obj.assemble_Kop_(obj.scatterer.dom_vie, freq);
            end
            
            % generate N operator for Ext domain
            if ~isequal(obj.setget_PrevDomain_ext(), obj.scatterer.dom_ext.r) || ...
                isempty(obj.operators.N_ext_) || (obj.setgetPrevFreq_VIE() ~= freq)
                
                % assemble operators for the Ext domain
                obj.operators.N_ext_ = obj.assemble_Nop_(obj.scatterer.dom_ext, freq);
                
                % save the operator to persistent variable
                obj.setget_Prev_N_ext(obj.operators.N_ext_);
                obj.setget_PrevDomain_ext(obj.scatterer.dom_ext.r);
                obj.setgetPrevFreq_VIE(freq);
                                
            end
            
            obj.operators.K_ext_ = obj.assemble_Kop_(obj.scatterer.dom_ext, freq);

            % assemble operators for the Near domain
            obj.operators.N_near_ = obj.assemble_Nop_(obj.scatterer.dom_near, freq);
            obj.operators.K_near_ = obj.assemble_Kop_(obj.scatterer.dom_near, freq);

            % get operator dimensions
            obj.get_operator_dimensions_vsie_();
                        
            % assemble projection matrices  
            obj.assemble_projection_matrix_(freq);

            % if gpu flag is nonzero allocate gpu device
            if gpu_flag
                
                % allocate operators N operators for vie and extended domain on GPU
                obj.operators.N_ext_gpu  = src_utils.to_GPU(obj.operators.N_ext_, tucker_flag);
                obj.operators.K_ext_gpu  = src_utils.to_GPU(obj.operators.K_ext_, tucker_flag);
                obj.operators.N_near_gpu = src_utils.to_GPU(obj.operators.N_near_, tucker_flag);
                obj.operators.K_near_gpu = src_utils.to_GPU(obj.operators.K_near_, tucker_flag);
                % move material properties to GPU
                obj.scatterer.prop_vie.Mcr = src_utils.to_GPU(obj.scatterer.prop_vie.Mcr, 0);
                obj.scatterer.prop_vie.Mrc = src_utils.to_GPU(obj.scatterer.prop_vie.Mrc, 0);
                    
            end
            
        end
    
        % ---------------------------------------------------------------------- %
        
        function assemble_final_mvp_(obj, freq)
            % private method assemble_final_mvp_() sets up the final VSIE
            % mvp for iterative solver
            
            % call mvp method to set final vsie mvp
            obj.mvp.set_final_mvp_vsie(obj.precorrection, obj.coil, freq);
        end
        
        % ---------------------------------------------------------------------- %
        
        function form_rhs(obj, freq)
            % method form_rhs() 
            
            solver_mode   = obj.task_settings_.vsie.Solver_mode;
            coupling_mode = obj.task_settings_.vsie.Coupling_Mode;
            
            
            if strcmp(solver_mode, 'Coupled')
                
                % allocate memory for rhs
                obj.rhs = zeros(obj.dims.N_sie + obj.dims.N_scat * obj.dims.ql, obj.dims.N_feeds);
                
                rhs_SIE = obj.rhs_c;
                
                % add excitation due to applied voltage
                obj.rhs(1:obj.dims.N_sie,:) = rhs_SIE;
                
            else
                
                rhs_SIE = obj.rhs_c;
                
                %% ????
                obj.rhs_c = obj.Zcoil_inv_ * rhs_SIE;
                
                Jc_ini = obj.operators.Jc_ini;
                
                switch coupling_mode
                    
                    case 'Basis'
                        
                        obj.rhs_b = gather(obj.operators.M_q2ql * (obj.operators.U_N * (obj.operators.X_N * rhs_SIE)));
                        
                    case 'Dense'
                        
                        obj.rhs_b = gather(obj.operators.Zbc_Nop * Jc_ini);
                        
                    case 'pFFT'
                        
                        gpu_flag = obj.task_settings_.general.GPU_flag;
                        emu = src_utils.EM_utils(freq);
                        
                        U = obj.precorrection.Z_C2B.Nop  * Jc_ini;
                        S = obj.projector.PS(:,obj.dims.N_sie+1:end);
                        P = obj.projector.PS(:,1:obj.dims.N_sie);
                        
                        if gpu_flag
                            Jin = gpuArray(P * Jc_ini);
                            L_ext = @(Jin) obj.mvp.L_ext_gpu(Jin);
                        else
                            Jin = P * Jc_ini;
                            L_ext = @(Jin) obj.mvp.L_ext(Jin);
                        end
                        
                        for i =1:size(Jin,2)
                            
                            E = L_ext(Jin(:,i));
                            Q(:,i) = 1 / emu.ce * S.' * E(:);
                        end
                        
                        obj.rhs_b = gather(U + Q);
                end
            end
        end
        
        
        
        % ---------------------------------------------------------------------- %
        
        function get_operator_dimensions_vsie_(obj)
            
            % Get vsie-related operator dimensions
            obj.get_operator_dimensions_vie_();
            
            % get operator dimensions
            obj.dims.op_ext  = obj.get_operator_dimensions_(obj.operators.N_ext_);
            obj.dims.op_near = obj.get_operator_dimensions_(obj.operators.N_near_);
        end        

    
        % ---------------------------------------------------------------------- %
        
        function generate_basis_(obj, freq)
           % checks if basis exists 
           
           % get path and filename 
           path     = obj.task_settings_.basis.Path;
           filename = obj.task_settings_.basis.Filename; 
           
           file = fullfile(path, filename);
           
           if 2 ~= exist(file,'file')
               src_assembly.generate_basis(obj.task_settings_, obj.dims, ...
                                           obj.scatterer, obj.coil, freq);
           end
           
        end
        
        
        % ---------------------------------------------------------------------- %
        function assemble_basis_coupling_(obj, freq)
            
            coup_type = obj.task_settings_.basis.Type;

            switch coup_type
                
                case 'Cross'
                    obj.coupling = Coupling_Basis_Cross(obj.task_settings_, obj.scatterer, obj.coil, obj.operators,...
                                                        obj.dims, obj.Zcoil_inv_, freq);
                case 'Dense'
                    obj.coupling = Coupling_Basis_Dense(obj.task_settings_, obj.scatterer, obj.coil, obj.operators,...
                                                        obj.dims, obj.Zcoil_inv_, freq);                             
                otherwise
                    error ('Error! Unknown mode of basis coupling!');
            end
            
            obj.coupling.assemble_basis_coupling();
            
        end

    end
    
% ====================================================================== %

end
