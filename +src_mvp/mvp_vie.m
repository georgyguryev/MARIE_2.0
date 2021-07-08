classdef mvp_vie < src_mvp.mvp_base
% Class mvp_vie sets up the VIE-related matrix-vector products 
% Author: Georgy Guryev, Cambridge, MA, 2019     
    
    % ====================================================================== %
    properties 
    
        PS
        precorrection   % precorrection object
    
    end
    % ====================================================================== %
    
    methods
        
        % define 
         function obj = mvp_vie(task_settings, operators, scatterer, projector, Z_coil, dims)
            
             % constructor for class mvp
            obj = obj@src_mvp.mvp_base(task_settings, operators, scatterer, dims);
                        
            % define core mvps (mvp.N, mvp.K, etc.)
            obj.set_core_mvps_vie_(projector);
            
            % set Gram mvp for preconditioner
            obj.set_precond_G_mvp_()

            % get impedance matrix from assembly
            obj.Z_coil = Z_coil;
                        
         end
         
         
         % ------------------------------------------------------------- %
         
         function set_final_mvp_vie(obj, precorrection, freq)
             % method set_vie_solver_mvp_() sets up the mvp for the vie
             % solver
            
            % get gpu flag
            gpu_flag = obj.task_settings_.general.GPU_flag;
            excitation_type = obj.task_settings_.vie.Excitation_type;
            rhs_filename    = obj.task_settings_.vie.Excitation_rhs;
            coupling_mode   = obj.task_settings_.vsie.Coupling_Mode;
            
            assemble_rhs_flag = (strcmp(excitation_type, 'Coil')) && ...
                                (~exist(rhs_filename));
                            
            S_ql = obj.scatterer.index_vie.index_ql(obj.dims.ql, obj.dims.Nvox_vie);
            
            if assemble_rhs_flag
                % setup coupling mvps
                switch coupling_mode
                    case 'pFFT'
                        
                        P  = obj.PS(:,1:obj.dims.N_sie);
                        S  = obj.PS(:,obj.dims.N_sie+1:end);
                        Jb = zeros(obj.dims.ql * obj.dims.N_scat, 1);
                        obj.c2b_coupling_vie_Nop_ = @(Jin) coupling_c2b_mvp(Jin, obj, precorrection, obj.dims, freq);
                        obj.c2b_coupling_vie_Kop_ = @(Jin) Jcb2Htot(Jin, Jb, @(Jtot)obj.K_ext(Jtot), precorrection.Z_C2B.Kop, P, S);
                    otherwise
                        obj.c2b_coupling_vie_Nop_ = @(Jin) coupling_c2b_dense_mvp(Jin, obj.operators_.Zbc_Nop);
                        obj.c2b_coupling_vie_Kop_ = @(Jin) coupling_c2b_dense_mvp(Jin, obj.operators_.Zbc_Kop);
                end
            end
            
            obj.Jb2Eb_mvp_     = @(Jin) Jb2Eb_mvp(Jin, obj.dims, obj.scatterer, freq);
            %             obj.implicit_b2c_coupling_ = @(Jin) src_mvp.coupling_b2c_mvp(Jin, obj, obj.precorrection, obj.dims, freq);
%             
%             obj.Jc2K_coupling_ = @(Jin) src_mvp.Jc2K_coupling(Jin,obj.scatterer.dom_vie.r, obj.precorrection.coil, obj.dims.ql, freq);
%             obj.Jb2Escat_mvp_  = @(Jin) src_mvp.Jb2Escat_mvp(Jin, obj, obj.dims, freq);

            % setup "VIE" mvp
            if gpu_flag
                obj.VIE_mvp_ = @(Jin) VIE_mvp_gpu(Jin, obj, obj.scatterer, obj.dims, S_ql, freq);                
            else
                obj.VIE_mvp_ = @(Jin) VIE_mvp(Jin, obj, obj.scatterer, obj.dims, S_ql, freq);
            end
          
         end


        
    end
    
    
  % ====================================================================== %
    
    methods (Access = private)
        
        % ------------------------------------------------------------- %
        
        function set_core_mvps_vie_(obj, projector)
            % function set_N_mvp_(obj) sets appropriate  mvps for vie
            % problem
            
%             obj.PS   = projector.PS;
            gpu_flag = obj.task_settings_.general.GPU_flag;
            excitation_type = obj.task_settings_.vie.Excitation_type;
            rhs_filename    = obj.task_settings_.vie.Excitation_rhs;
            coupling_mode   = obj.task_settings_.vsie.Coupling_Mode;
                                                 
                                                     
            [N_mvp_vie, K_mvp_vie, G_mvp_vie] = obj.set_core_mvps_(obj.operators_.N_vie_,...
                                                                   obj.operators_.K_vie_,...
                                                                   obj.dims.vie, obj.dims.op_vie);  
                                                                                                                                
            % save mvps to obj.mvp
            obj.N_mvp_vie_  = N_mvp_vie;
            obj.K_mvp_vie_  = K_mvp_vie;
            obj.G_mvp_vie_  = G_mvp_vie;
            
            if strcmp(coupling_mode, 'pFFT')
                
                [N_mvp_ext, K_mvp_ext, G_mvp_ext] = obj.set_core_mvps_(obj.operators_.N_ext_,...
                                                                       obj.operators_.K_ext_,...
                                                                       obj.dims.ext, obj.dims.op_ext);
                                                                   
                 % setup core mvps for near domain                                        
                [N_mvp_near, K_mvp_near, G_mvp_near] = obj.set_core_mvps_(obj.operators_.N_near_,...
                                                             obj.operators_.K_near_,...
                                                             obj.dims.near,...
                                                             obj.dims.op_near);                                                   
                                                                   
                obj.N_mvp_ext_  = N_mvp_ext;
                obj.K_mvp_ext_  = K_mvp_ext;
                obj.G_mvp_ext_  = G_mvp_ext;
                
                
                % for extended domain
                obj.N_mvp_near_ = N_mvp_near;
                obj.K_mvp_near_ = K_mvp_near;
                obj.G_mvp_near_ = G_mvp_near;
                
                % copy projection/interpolation matrices
                obj.PS = projector.PS;
                
                % setup L operators
                obj.L_mvp_ext_    = @(Jin) L_mvp_pwx(Jin, N_mvp_ext, G_mvp_ext, obj.dims.ext);
                obj.L_mvp_near_   = @(Jin) L_mvp_pwx(Jin, N_mvp_near, G_mvp_near, obj.dims.near);
                
                if gpu_flag
                    
                    % setup core mvps for extended domain
                    [N_mvp_ext_gpu, K_mvp_ext_gpu, G_mvp_ext_gpu] = obj.set_core_mvps_(obj.operators_.N_ext_gpu,...
                        obj.operators_.K_ext_gpu,...
                        obj.dims.ext, obj.dims.op_ext);
                    
                    
                    % setup core mvps for near domain
                    [N_mvp_near_gpu, K_mvp_near_gpu, G_mvp_near_gpu] = obj.set_core_mvps_(obj.operators_.N_near_gpu,...
                        obj.operators_.K_near_gpu,...
                        obj.dims.near,...
                        obj.dims.op_near);
                    
                    
                    obj.K_mvp_ext_gpu  = K_mvp_ext_gpu;
                    obj.K_mvp_near_gpu = K_mvp_near_gpu;
                    
                    % setup L operators on gpu
                    obj.L_mvp_ext_gpu  = @(Jin) L_mvp_pwx_gpu(Jin, N_mvp_ext_gpu, G_mvp_ext_gpu, obj.dims.ext);
                    obj.L_mvp_near_gpu = @(Jin) L_mvp_pwx_gpu(Jin, N_mvp_near_gpu, G_mvp_near_gpu, obj.dims.near);
                    
                    % setup projection and interpolation mvps on gpu
                    obj.Jtot2Jcb_mvp_gpu = @(Jtot) Jtot2Jcb_mvp(Jtot, gpuArray(obj.PS));
                    obj.Jcb2Jtot_mvp_gpu = @(Jcb)  Jcb2Jtot_mvp(Jcb, gpuArray(obj.PS));
                    
                end
               
            end
            
            % setup L operators
            obj.L_mvp_vie_    = @(Jin) L_mvp_pwx(Jin, N_mvp_vie, G_mvp_vie, obj.dims.vie);
        
            % setup projection and interpolation mvps
            obj.Jtot2Jcb_mvp_ = @(Jtot) Jtot2Jcb_mvp(Jtot, obj.PS);
            obj.Jcb2Jtot_mvp_ = @(Jcb)  Jcb2Jtot_mvp(Jcb, obj.PS);
                         
            % decoupled solver
            if gpu_flag
                             
                % setup core mvps for scatterer domain
                [N_mvp_vie_gpu, K_mvp_vie_gpu, ~] = obj.set_core_mvps_(obj.operators_.N_vie_gpu,...
                                                                       obj.operators_.K_vie_gpu,...
                                                                       obj.dims.vie, obj.dims.op_vie);  
                                                       
                obj.N_mvp_vie_gpu  = N_mvp_vie_gpu;
                obj.K_mvp_vie_gpu  = K_mvp_vie_gpu;
               
            end

        end
        
                                 
    end
    
   % ====================================================================== %
   
end