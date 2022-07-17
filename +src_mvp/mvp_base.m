classdef (Abstract) mvp_base < handle
% Abstract class mvp_base defines interface and shared properties for mvps  
% Author: Georgy Guryev, Cambridge, MA, 2019     
    
    properties %(Access = protected)
        
        % task settings
        task_settings_
        
        % N,K operators for various domains
        operators_

        % scatterer
        scatterer
        
        % coil impedance matrix
        Z_coil
        Z_coil_inv
        coil
        
        % hadles for mvp functions; set up in constructor
        G_mvp_vie_
        N_mvp_vie_
        N_mvp_vie_gpu
        K_mvp_vie_
        K_mvp_vie_gpu
        L_mvp_vie_
        L_mvp_vie_gpu
        
        G_mvp_ext_
        G_mvp_precond_
        N_mvp_ext_
        L_mvp_ext_
        L_mvp_ext_gpu
        K_mvp_ext_
        K_mvp_ext_gpu
        
        G_mvp_near_
        N_mvp_near_
        K_mvp_near_
        L_mvp_near_
        K_mvp_near_gpu
        L_mvp_near_gpu
        
        SIE_mvp_
        
        VIE_mvp_
        VIE_mvp_gpu
        Jb2Escat_mvp_
        Jb2Eb_mvp_
        Jtot2Jcb_mvp_
        Jcb2Jtot_mvp_
        Jtot2Jcb_mvp_gpu
        Jcb2Jtot_mvp_gpu
        
        implicit_c2b_coupling_
        implicit_b2c_coupling_
        implicit_UWU_coupling_
        
        c2b_coupling_vie_
        c2b_coupling_vie_Nop_
        c2b_coupling_vie_Kop_
        
        
        Jc2K_coupling_ 
        Jc2H_coupling_
        Jb2H_coupling_
        Jb2H_tot_
        
        VSIE_mvp_ % handle to final matrix-vector product
        VSIE_mvp_temp_

                
        dims
        
    end
    
    % ====================================================================== %
 
 
    methods 
        function obj = mvp_base(task_settings, operators, scatterer, dims)
            % constructor of abstract class mvp_base
            obj.task_settings_ = task_settings;
            obj.operators_     = operators;
            obj.scatterer     = scatterer;
            obj.dims = dims;
        end        
        
               
        % ------------------------------------------------------------- %
      
        function Jout = G_vie(obj,Jin)
            Jout = obj.G_mvp_vie_(Jin);
        end
        
        % ------------------------------------------------------------- %

        function Jout = N_vie(obj,Jin)
            Jout = obj.N_mvp_vie_(Jin);
        end
        
        % ------------------------------------------------------------- %
        
        function Jout = N_vie_gpu(obj,Jin)
            Jout = obj.N_mvp_vie_gpu(Jin);
        end
        
         % ------------------------------------------------------------- %

        function Jout = K_near(obj,Jin)
            Jout = obj.K_mvp_near_(Jin);
        end
        
        % ------------------------------------------------------------- %

        function Jout = K_near_gpu(obj,Jin)
            Jout = obj.K_mvp_near_gpu(Jin);
        end
        
        % ------------------------------------------------------------- %

        function Jout = K_vie(obj,Jin)
            Jout = obj.K_mvp_vie_(Jin);
        end
        
        
        % ------------------------------------------------------------- %

        function Jout = K_vie_gpu(obj,Jin)
            Jout = obj.K_mvp_vie_gpu(Jin);
        end
        % ------------------------------------------------------------- %

        function Jout = K_ext(obj,Jin)
            Jout = obj.K_mvp_ext_(Jin);
        end
        
        
        % ------------------------------------------------------------- %

        function Jout = K_ext_gpu(obj,Jin)
            Jout = obj.K_mvp_ext_gpu(Jin);
        end
        % ------------------------------------------------------------- %

        function Jout = L_vie(obj,Jin)
            Jout = obj.L_mvp_vie_(Jin);
        end
        
        
        % ------------------------------------------------------------- %

        function Jout = L_vie_gpu(obj,Jin)
            Jout = obj.L_mvp_vie_gpu(Jin);
        end
        
        
        % ------------------------------------------------------------- %
        function Jout = L_ext(obj,Jin)
            Jout = obj.L_mvp_ext_(Jin);
        end
        
        % ------------------------------------------------------------- %

        function Jout = L_ext_gpu(obj,Jin)
            Jout = obj.L_mvp_ext_gpu(Jin);
        end
        
        % ------------------------------------------------------------- %
        function Jout = L_near(obj, Jin)
            Jout = obj.L_mvp_near_(Jin);
        end
        
        % ------------------------------------------------------------- %

        function Jout = L_near_gpu(obj, Jin)
            Jout = obj.L_mvp_near_gpu(Jin);
        end
        % ------------------------------------------------------------- %
            
        function Jout = Jb2Escat(obj,Jin)
            Jout = obj.Jb2Escat_mvp_(Jin);
        end
       % ------------------------------------------------------------- %
            
        function Jout = Jb2Eb(obj,Jin)
            Jout = obj.Jb2Eb_mvp_(Jin);
        end
        % ------------------------------------------------------------- %
        
        function Jout = VIE(obj, Jin)
            Jout = obj.VIE_mvp_(Jin);
        end
        
        % ------------------------------------------------------------- %
        
        function [Icu, Ics] = SIE(obj, Jbu, Jbs, Z_L, rhs_cp)
            [Icu, Ics] = obj.SIE_mvp_(Jbu, Jbs, Z_L, rhs_cp);
        end
        % ------------------------------------------------------------- %
        
        function Jout = VIE_gpu(obj, Jin)
            Jout = obj.VIE_mvp_gpu(Jin);
        end
        
        % ------------------------------------------------------------- %

        function Jout = VSIE(obj, Jin)
                Jout = obj.VSIE_mvp_(Jin);
        end
        
        % ------------------------------------------------------------- %
        
        function Jout = Jtot2Jcb(obj, Jin)
            Jout = obj.Jtot2Jcb_mvp_(Jin);
        end
        
        % ------------------------------------------------------------- %

        function Jout = Jcb2Jtot(obj,Jin)
            Jout = obj.Jcb2Jtot_mvp_(Jin);
        end
        
        % ------------------------------------------------------------- %
        
        function Jout = Jtot2Jcb_gpu(obj, Jin)
            Jout = obj.Jtot2Jcb_mvp_gpu(Jin);
        end
        
        % ------------------------------------------------------------- %

        function Jout = Jcb2Jtot_gpu(obj,Jin)
            Jout = obj.Jcb2Jtot_mvp_gpu(Jin);
        end
        
        % ------------------------------------------------------------- %
        
        function Jout = G_prec(obj, Jin)
            Jout = obj.G_mvp_precond_(Jin);
        end
                
        % ------------------------------------------------------------- %
        
        function Jout = c2b_implicit_coupling(obj,Jin)
            Jout = obj.implicit_c2b_coupling_(Jin);
        end
        
        
        function Einc = c2b_coupling_vie_Nop(obj, Jin)
            Einc = obj.c2b_coupling_vie_Nop_(Jin);
        end
        
        % ------------------------------------------------------------- %
        
        function Hinc = c2b_coupling_vie_Kop(obj, Jin)
            Hinc = obj.c2b_coupling_vie_Kop_(Jin);
        end
        
        % ------------------------------------------------------------- %
        
        function Jout = b2c_implicit_coupling(obj, Jin)
            Jout = obj.implicit_b2c_coupling_(Jin);
        end
        
        % ------------------------------------------------------------- %
        
        function Jout = uwu_implicit_coupling(obj, Jin)
            Jout = obj.implicit_UWU_coupling_(Jin);
        end
        
        % ------------------------------------------------------------- %
        
        function Jout = get_Hinc_coupling(obj, Jin, Ic0, Y_LC, rhs_cp, type)
            
            switch type
                case 'body'
                    Jout = obj.Jb2H_coupling_(Jin, Ic0, Y_LC, rhs_cp);

                case 'coil'
                    Jout = obj.Jc2H_coupling_(Jin);
            end
        end 
        
        % ------------------------------------------------------------- %
        
        function Htot = Jb2Htot(obj, Jb)
            Htot = obj.Jb2H_tot_(Jb);
        end

    end
    
  % ====================================================================== %
  methods (Access = protected)
      
      
      function [N_mvp, K_mvp, G_mvp] = set_core_mvps_(obj, Nop, Kop, task_dims, op_dims)
          % private method set_mvps_core_() is used to set fundamental mvps for the
          % vie/vsie task.

          
          % get flags with setting
          flag_pwx    = obj.task_settings_.vie.PWX;
          flag_tucker = obj.task_settings_.vie.Tucker;
          
          % set appropriate implementation of mvp product
          switch (2 * flag_pwx + flag_tucker)
              
              case 0
                  N_mvp = @(Jin) N_mvp_pwc(Jin, Nop, task_dims, op_dims);
                  K_mvp = @(Jin) K_mvp_pwc(Jin, Kop, task_dims, op_dims);
                  G_mvp = @(Jin) G_mvp_pwc(Jin, obj.scatterer.dom_vie.res, task_dims);
              case 1
                  N_mvp = @(Jin) N_mvp_pwc_tucker(Jin, Nop, task_dims, op_dims);
                  K_mvp = @(Jin) K_mvp_pwc_tucker(Jin, Kop, task_dims, op_dims);
                  G_mvp = @(Jin) G_mvp_pwc(Jin, obj.scatterer.dom_vie.res, task_dims);
              case 2
                  N_mvp = @(Jin) N_mvp_pwl(Jin, Nop, task_dims, op_dims, obj.dims);
                  K_mvp = @(Jin) K_mvp_pwl(Jin, Kop, task_dims, op_dims);
                  G_mvp = @(Jin) G_mvp_pwl(Jin,obj.scatterer.dom_vie.res, task_dims, obj.dims);
                  
              case 3
                  N_mvp = @(Jin) N_mvp_pwl_tucker(Jin, Nop, task_dims, op_dims);
                  K_mvp = @(Jin) K_mvp_pwl_tucker(Jin, Kop, task_dims, op_dims);
                  G_mvp = @(Jin) G_mvp_pwl(Jin,obj.scatterer.dom_vie.res, task_dims, obj.dims);
              otherwise
                  error('Wrong operator format!');
          end
      end
      
      function set_precond_G_mvp_(obj)
          % method set_precond_G_mvp_() sets up the Grammian
          % multiplication for preconditioner (unlike regular one, it is
          % aplied once at the assembly stage to a vector of material
          % properties for body voxels
          obj.G_mvp_precond_ = @(Jin) G_mvp_precond(Jin, obj.scatterer.dom_vie.res, obj.dims);
      end
      
      
  end
    
end