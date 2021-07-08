classdef mvp_vsie < src_mvp.mvp_base
% Class mvp_vsie sets up the VSIE-related matrix-vector products 
% Author: Georgy Guryev, Cambridge, MA, 2019     
    
% ====================================================================== %
   
    properties (Access = private)
                
        PS              % projection matrix
        precorrection   % precorrection object
    end
    
% ====================================================================== %
    
    methods
        
        % define
        function obj = mvp_vsie(task_settings, operators, scatterer, projector, Z_coil, Z_coil_inv, dims)
            
            % constructor for class mvp
            obj = obj@src_mvp.mvp_base(task_settings, operators, scatterer, dims);
            
            % define mvp_N and mvp_K
            obj.set_core_mvps_vsie_(projector,scatterer.freq);
            
            % define G_mvp_precond
            obj.set_precond_G_mvp_()
            
            % get impedance matrix from assembly
            obj.Z_coil = Z_coil;
            obj.Z_coil_inv = Z_coil_inv;
        end
         
        % ------------------------------------------------------------- %
         
        function set_final_mvp_vsie(obj, precorrection, coil, freq)
            % method set_final_mvp_vsie(obj) sets up the final mvp for the
            % iterative solver
            
            obj.coil = coil;

            % get solver and coupling modes
            coupling_type = obj.task_settings_.vsie.Coupling_Mode;   

            switch coupling_type
                
                case 'pFFT'
                    % copy precorrection matrices
                    obj.precorrection = precorrection;
                    
                    % assemble core mvps
                    obj.set_final_mvp_pFFT_(coil, freq);
                    
                case 'Basis'
                    obj.set_final_mvp_basis_(coil, freq);
                    
                case 'Dense'
                    % assemble core mvps
                    obj.set_final_mvp_explicit_(coil, freq);
            end
            
        end
        
        % ------------------------------------------------------------- %
        
        function update_final_mvp_decoupled_implicit(obj, vie_solver)
            % update_final_mvp_decoupled_implicit(vie_solver)
            
            obj.VSIE_mvp_ = @(Jcb)obj.VSIE_mvp_temp_(Jcb, vie_solver);
        end
        
    end
    
  % ====================================================================== %
    
    methods (Access = private)
        
        % ------------------------------------------------------------- %
        function set_final_mvp_pFFT_(obj, coil, freq)
            
            % get gpu flag
            gpu_flag    = obj.task_settings_.general.GPU_flag;
            Solver_mode = obj.task_settings_.vsie.Solver_mode;

            % setup coupling mvps
            obj.implicit_c2b_coupling_ = @(Jin) coupling_c2b_mvp(Jin, obj, obj.precorrection, obj.dims, freq);
            
            % 
            obj.Jc2K_coupling_ = @(Jin) Jc2K_coupling(Jin,obj.scatterer.dom_vie.r, coil, obj.dims.ql, freq);
            obj.Jc2H_coupling_ = @(Jin) Jc2H_coupling_basis(Jin, obj.operators_.U_K, obj.operators_.alpha_K,...
                                                                    obj.operators_.M_q2ql);
            obj.Jb2Escat_mvp_  = @(Jin) Jb2Escat_mvp(Jin, obj, obj.dims, freq);
            obj.Jb2Eb_mvp_     = @(Jin) Jb2Eb_mvp(Jin, obj.dims, obj.scatterer, freq);

            % get problem sizes 
            N_b   = obj.dims.N_scat * obj.dims.ql;
            N_sie = obj.dims.N_sie;

            % get electro-magnetic parameters  
            emu = src_utils.EM_utils(freq);
            
            
            if strcmp(Solver_mode, 'Coupled')
            
                % define mixing matrices (for optimized of the final mvp)
                P1 = 1 / emu.ce * blkdiag(speye(N_sie, N_sie), -speye(N_b, N_b));
                P2 = 1 / emu.ce * blkdiag(sparse(N_sie, N_sie), speye(N_b, N_b));
                P3 = speye(N_sie + N_b);

                % resulting mixing matrix
                P_tot = [P1, P2, P3];

                if gpu_flag

                    P_tot  = gpuArray(P_tot);
                    obj.PS = gpuArray(obj.PS);
                    obj.precorrection.Z_tot = gpuArray(obj.precorrection.Z_tot);

                    obj.VSIE_mvp_ = @(Jcb) VSIE_coupled_pFFT_gpu(Jcb, obj,...
                                                                     obj.scatterer,...
                                                                     obj.dims,...
                                                                     obj.precorrection,...
                                                                     P_tot);
                else

                    % set mvp for sie
                    obj.VSIE_mvp_ = @(Jcb) VSIE_coupled_pFFT(Jcb, obj,...
                                                                     obj.scatterer,...
                                                                     obj.dims,...
                                                                     obj.precorrection,...
                                                                     P_tot);           
                end
                
            else
                
                P = obj.PS(:,1:obj.dims.N_sie);
                S = obj.PS(:,obj.dims.N_sie+1:end);
                
                if gpu_flag
                    P = src_utils.to_GPU(P,0);
                    S = src_utils.to_GPU(S,0);
                    Z_C2B = src_utils.to_GPU(obj.precorrection.Z_C2B.Nop,0);
                    Z_C2B_Kop = src_utils.to_GPU(obj.precorrection.Z_C2B.Kop,0);
                    Z_C2B_T = Z_C2B.';
                    obj.Z_coil_inv = src_utils.to_GPU(obj.Z_coil_inv,0);
                    Mc_inv_Gram = gpuArray(obj.scatterer.prop_ext.Mc_inv * obj.scatterer.dom_ext.res.^3);
                    L_ext = @(Jin) obj.L_ext_gpu(Jin);
                    K_ext = @(Jin) obj.K_ext_gpu(Jin);
                else
                    Z_C2B     = obj.precorrection.Z_C2B.Nop;
                    Z_C2B_Kop = obj.precorrection.Z_C2B.Kop;
                    Z_C2B_T = Z_C2B.';
                    Mc_inv_Gram = obj.scatterer.prop_ext.Mc_inv * obj.scatterer.dom_ext.res.^3;
                    L_ext = @(Jin) obj.L_ext(Jin);
                    K_ext = @(Jin) obj.K_ext(Jin);
                end
                
                obj.VSIE_mvp_ = @(Jcb) VSIE_decoupled_pFFT(Jcb, L_ext, obj.dims, Mc_inv_Gram,...
                                                                 Z_C2B, Z_C2B_T, obj.Z_coil_inv,...
                                                                 P, S, freq);
                                                             
                obj.Jcb2H_tot_ = @(Jc, Jb) Jcb2Htot(Jc, Jb, K_ext, Z_C2B_Kop, P,S);
                                                             
                obj.SIE_mvp_ = @(Jb, port) SIE_decoupled_pFFT(Jb, port, obj.operators_.Jc_ini,...
                                                              obj.Z_coil_inv, L_ext, P, S, freq);                                             
                                                             
            end
        end
               
        % ------------------------------------------------------------- %
        
        function set_final_mvp_basis_(obj, coil, freq)
            % function set_final_mvp_coupled_implicit_();
            
            Solver_mode = obj.task_settings_.vsie.Solver_mode;
            % 
            obj.Jc2K_coupling_ = @(Jin) Jc2K_coupling(Jin,obj.scatterer.dom_vie.r, coil, obj.dims.ql, freq);
            obj.Jc2H_coupling_ = @(Jin) Jc2H_coupling_basis(Jin, obj.operators_.U_K, obj.operators_.alpha_K,...
                                                                    obj.operators_.M_q2ql);
            obj.Jb2Escat_mvp_  = @(Jin) Jb2Eb_vie_mvp(Jin, obj, freq);
            obj.Jb2Eb_mvp_     = @(Jin) Jb2Eb_mvp(Jin, obj.dims, obj.scatterer, freq);

            if strcmp(Solver_mode, 'Coupled')
                obj.VSIE_mvp_ = @(Jcb) VSIE_coupled_basis(Jcb, obj, obj.dims);
            else
                obj.VSIE_mvp_ = @(Jb) VSIE_decoupled_basis(Jb, obj);
                obj.SIE_mvp_  = @(Jb, port) SIE_decoupled_basis(Jb, port, obj.operators_.U_N,...
                                                                        obj.operators_.X_N,...
                                                                        obj.operators_.Jc_ini,...
                                                                        obj.operators_.M_q2ql);
            end
        end
        
        % ------------------------------------------------------------- %
        function set_final_mvp_explicit_(obj, coil, freq)
            
            Solver_mode = obj.task_settings_.vsie.Solver_mode;
            gpu_flag_coupling = obj.task_settings_.basis.GPU_flag;
            
            obj.Jc2H_coupling_ = @(Jin) Jc2H_coupling_dense(Jin, obj.operators_.Zbc_Kop);
            obj.Jb2Escat_mvp_  = @(Jin) Jb2Eb_vie_mvp(Jin, obj, freq);
            obj.Jb2Eb_mvp_     = @(Jin) Jb2Eb_mvp(Jin, obj.dims, obj.scatterer, freq);

            if strcmp(Solver_mode, 'Coupled')

                obj.VSIE_mvp_ = @(Jcb) VSIE_coupled_explicit(Jcb, obj, obj.Z_coil,...
                                                                     obj.operators_.Zbc_Nop, obj.operators_.Zbc_Nop.');
            else 
                % set up decoupled vsie mvp (solve for tissue currents)
%                 obj.VSIE_mvp_ = @(Jcb) src_mvp.VSIE_decoupled_explicit(Jcb, obj, obj.Z_coil,...
%                                                                        obj.operators_.Zbc, obj.operators_.Zbc.');
                if gpu_flag_coupling
                    obj.Z_coil_inv = gpuArray(obj.Z_coil_inv);
                end
                
                obj.VSIE_mvp_ = @(Jcb) VSIE_decoupled_explicit(Jcb, obj, obj.Z_coil_inv,...
                                                                       obj.operators_.Zbc_Nop);
                                                                   
                % set up decoupled sie mvp                                                
                obj.SIE_mvp_ = @(Jb, port) SIE_decoupled_explicit(Jb, obj.operators_.Jc_ini,...
                                                                          obj.Z_coil,...
                                                                          obj.operators_.Zbc_Nop,...
                                                                          port);                                               
            end                                            

        end
        % ------------------------------------------------------------- %

        
        function set_core_mvps_vsie_(obj, projector, freq)
            % function set_N_mvp_(obj) sets appropriate  mvps for vie
            % problem
            
            
            % get solver and coupling modes
            coupling_type = obj.task_settings_.vsie.Coupling_Mode;
            
            % select mvp assembly type based on assembly_type
            switch coupling_type
                
                % coupled implicit
                case 'Dense'
                    
                    obj.set_core_mvps_explicit_(freq);
                    
                % coupled explicit
                case 'pFFT'
                    
                     % copy projection/interpolation matrices
                    obj.PS = projector.PS;
                    
                    % assemble core mvps 
                    obj.set_core_mvps_pFFT_();
                    
                % decoupled implicit    
                case 'Basis' 
                    
                    % assemble core mvps
                    obj.set_core_mvps_basis_(freq);
                    
            end

        end
        
        % ------------------------------------------------------------- %
        
        function set_near_mvp_(obj)
            % method set_near_mvp_() sets up the neat mvp for mvp
            obj.N_near_mvp = @(Jin) N_mvp_pwc(Jin, obj.N_mvp_near_, obj.dims.near, op_dims);
        end
        
        
        % ------------------------------------------------------------- %
       
        
        function set_core_mvps_pFFT_(obj)
            % function set_core_mvps_coupled_implicit_(obj)
            
            gpu_flag = obj.task_settings_.general.GPU_flag;
            
            % setup core mvps for extended domain
            [N_mvp_ext, K_mvp_ext, G_mvp_ext] = obj.set_core_mvps_(obj.operators_.N_ext_,...
                                                         obj.operators_.K_ext_,...
                                                         obj.dims.ext, obj.dims.op_ext);
                                                     

            [N_mvp_vie, K_mvp_vie, G_mvp_vie] = obj.set_core_mvps_(obj.operators_.N_vie_,...
                                                                   obj.operators_.K_vie_,...
                                                                   obj.dims.vie, obj.dims.op_vie);  
                                                               
                                                      
            % setup core mvps for near domain                                        
            [N_mvp_near, K_mvp_near, G_mvp_near] = obj.set_core_mvps_(obj.operators_.N_near_,...
                                                             obj.operators_.K_near_,...
                                                             obj.dims.near,...
                                                             obj.dims.op_near);
                                                                  
            % store mvp hadles in object properties
            obj.N_mvp_ext_  = N_mvp_ext;
            obj.G_mvp_ext_  = G_mvp_ext;
            obj.K_mvp_ext_  = K_mvp_ext;
            
            obj.N_mvp_vie_  = N_mvp_vie;
            obj.K_mvp_vie_  = K_mvp_vie;
            obj.G_mvp_vie_  = G_mvp_vie;
            
            % for extended domain          
            obj.N_mvp_near_ = N_mvp_near;
            obj.K_mvp_near_ = K_mvp_near; 
            obj.G_mvp_near_ = G_mvp_near;
            
            % setup L operators 
            obj.L_mvp_ext_    = @(Jin) L_mvp_pwx(Jin, N_mvp_ext, G_mvp_ext, obj.dims.ext);
            obj.L_mvp_near_   = @(Jin) L_mvp_pwx(Jin, N_mvp_near, G_mvp_near, obj.dims.near);
            
            % setup projection and interpolation mvps
            obj.Jtot2Jcb_mvp_ = @(Jtot) Jtot2Jcb_mvp(Jtot, obj.PS);
            obj.Jcb2Jtot_mvp_ = @(Jcb)  Jcb2Jtot_mvp(Jcb, obj.PS);
            
            if gpu_flag
                
                % setup core mvps for extended domain
                [N_mvp_ext_gpu, K_mvp_ext_gpu, G_mvp_ext_gpu] = obj.set_core_mvps_(obj.operators_.N_ext_gpu,...
                                                                                   obj.operators_.K_ext_gpu,...
                                                                                   obj.dims.ext, obj.dims.op_ext);
                 
                obj.K_mvp_ext_gpu = K_mvp_ext_gpu;
                
                [N_mvp_vie_gpu, K_mvp_vie_gpu , ~] = obj.set_core_mvps_(obj.operators_.N_vie_gpu,...
                                                                        obj.operators_.K_vie_gpu,...
                                                                        obj.dims.vie, obj.dims.op_vie);                                      
                % setup core mvps for near domain
                [N_mvp_near_gpu, K_mvp_near_gpu, G_mvp_near_gpu] = obj.set_core_mvps_(obj.operators_.N_near_gpu,...
                                                                        obj.operators_.K_near_gpu,...
                                                                        obj.dims.near,...
                                                                         obj.dims.op_near);
                                                             
                obj.N_mvp_vie_gpu  = N_mvp_vie_gpu;
                obj.K_mvp_vie_gpu  = K_mvp_vie_gpu;
                obj.K_mvp_near_gpu = K_mvp_near_gpu;
                
                % setup L operators on gpu
                obj.L_mvp_ext_gpu  = @(Jin) L_mvp_pwx_gpu(Jin, N_mvp_ext_gpu, G_mvp_ext_gpu, obj.dims.ext);
                obj.L_mvp_near_gpu = @(Jin) L_mvp_pwx_gpu(Jin, N_mvp_near_gpu, G_mvp_near_gpu, obj.dims.near);
                
                % setup projection and interpolation mvps on gpu
                obj.Jtot2Jcb_mvp_gpu = @(Jtot) Jtot2Jcb_mvp(Jtot, gpuArray(obj.PS));
                obj.Jcb2Jtot_mvp_gpu = @(Jcb)  Jcb2Jtot_mvp(Jcb, gpuArray(obj.PS));
                
            end
        end
        
        % ------------------------------------------------------------- %
        
        function set_core_mvps_explicit_(obj, freq)
            
            % get gpu flag
            gpu_flag = obj.task_settings_.general.GPU_flag;
            
            % setup core mvps for body domain
            [N_mvp, K_mvp, G_mvp] = obj.set_core_mvps_(obj.operators_.N_vie_,...
                                                       obj.operators_.K_vie_,...
                                                       obj.dims.vie, obj.dims.op_vie);

            S_ql = obj.scatterer.index_vie.index_ql(obj.dims.ql, obj.dims.Nvox_vie);
                                     
            
            % store mvp hadles in object properties
            obj.N_mvp_vie_  = N_mvp;
            obj.K_mvp_vie_  = K_mvp;
            obj.G_mvp_vie_  = G_mvp;
            
            % setup "VIE" mvp
            if gpu_flag
                
                [N_mvp_vie_gpu, K_mvp_vie_gpu , ~] = obj.set_core_mvps_(obj.operators_.N_vie_gpu,...
                                                                        obj.operators_.K_vie_gpu,...
                                                                        obj.dims.vie, obj.dims.op_vie);
                                                                    

                obj.VIE_mvp_ = @(Jin) VIE_mvp_gpu(Jin, obj, obj.scatterer, obj.dims, S_ql, freq);
                
                obj.N_mvp_vie_gpu = N_mvp_vie_gpu;
                obj.K_mvp_vie_gpu = K_mvp_vie_gpu; 
                
            else
                obj.VIE_mvp_ = @(Jin) VIE_mvp(Jin, obj, obj.scatterer, obj.dims, S_ql, freq);
            end
            

        end
        
        % ----------------------------------------------------------------- %
        function set_core_mvps_basis_(obj,freq)
            
            % get gpu flag
            gpu_flag = obj.task_settings_.general.GPU_flag;
            
            Volume = obj.scatterer.dom_vie.res.^3;
            S_ql = obj.scatterer.index_vie.index_ql(obj.dims.ql, obj.dims.Nvox_vie);

        
            % setup coupling mvps
            obj.implicit_c2b_coupling_ = @(Jin) coupling_c2b_basis_mvp(Jin, obj.operators_.U_N, obj.operators_.alpha_N, Volume, obj.operators_.M_q2ql);
            obj.implicit_b2c_coupling_ = @(Jin) coupling_b2c_basis_mvp(Jin, obj.operators_.U_N, obj.operators_.alpha_N, Volume, obj.operators_.M_q2ql);
%             obj.implicit_UWU_coupling_ = @(Jin) coupling_UWU_basis_mvp(Jin, obj.operators_.U_N, obj.operators_.W_N, obj.operators_.M_q2ql);  
            obj.implicit_UWU_coupling_ = @(Jin) coupling_UWU_basis_mvp(Jin, obj.operators_.UU_N, obj.operators_.V_N, obj.operators_.M_q2ql);  

            [N_mvp_vie, K_mvp_vie, G_mvp_vie] = obj.set_core_mvps_(obj.operators_.N_vie_,...
                                                                   obj.operators_.K_vie_,...
                                                                   obj.dims.vie, obj.dims.op_vie);
                                                               
                                                      
                                                               
            % save matrix-vector products
            obj.N_mvp_vie_  = N_mvp_vie;
            obj.K_mvp_vie_  = K_mvp_vie;
            obj.G_mvp_vie_  = G_mvp_vie;
            
            obj.VIE_mvp_ = @(Jin) VIE_mvp(Jin, obj, obj.scatterer, obj.dims, S_ql, freq);


            % setup "VIE" mvp
            if gpu_flag
                
                [N_mvp_vie_gpu, K_mvp_vie_gpu , ~] = obj.set_core_mvps_(obj.operators_.N_vie_gpu,...
                                                                        obj.operators_.K_vie_gpu,...
                                                                        obj.dims.vie, obj.dims.op_vie);
                                                                    
                obj.VIE_mvp_ = @(Jin) VIE_mvp_gpu(Jin, obj, obj.scatterer, obj.dims, S_ql, freq);                
                
                obj.N_mvp_vie_gpu = N_mvp_vie_gpu;
                obj.K_mvp_vie_gpu = K_mvp_vie_gpu; 
                
            end
            
            
        end
        
      
        
        
    end
    
  % ====================================================================== %

end