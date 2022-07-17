classdef Precorrection < handle
    % definition of class precorrection - a class for
    % precorrection of C2B and C2C interactions 

% ====================================================================== %

    properties
        
        Z_C2C   % coil to coil precorrection         
        Z_C2B   % coil to body precorrection
        Z_tot   % total precorrection: [Z_C2C, Z_C2B; -Z_C2B.', 0]
    end
    
% ====================================================================== %
    
    properties (SetAccess=immutable)
        
        % system parameters
        coil
        scatterer
        projector
        task_settings
        dims
        mvp
        freq
        
        % Coil2Coil impedance matrix
        Zc
    end
    
% ====================================================================== %


    methods
        
        function obj = Precorrection(coil, scatterer, projector,...
                                     task_settings, dims, mvp, Zc, freq)
                 % constructor of Precorrection class 
                 
                 % store references to input parameters
                 obj.coil      = coil;
                 obj.scatterer = scatterer;
                 obj.projector = projector;
                 obj.dims      = dims;
                 obj.mvp       = mvp;
                 obj.freq      = freq;
                 obj.Zc        = Zc;
                 obj.task_settings = task_settings;

        end
        
        % --------------------------------------------------------------- %
        
        function assemble_precorrection(obj)
            % function assemble_precorrection() assembles precorrection for
            % Coil2Body and Coil2Coil interactions
            
            N_b = obj.dims.N_scat * obj.dims.ql;
            
            % get previous coil model
            prev_coil = obj.setgetPrevCoil();
   
            
            % assemble Coil2Body precorrection 
            obj.assemble_C2B_precorrection_();
            
            % get Coil2Coil precorrection 
            obj.Z_C2C =  obj.setget_Prev_Zc_prec();
            
            if strcmp('Coupled', obj.task_settings.vsie.Solver_mode)
                
                if (isempty(prev_coil) || nnz(prev_coil.elem ~= obj.coil.elem))
                    
                    
%                     % assemble Coil2Coil precorrection
%                     obj.assemble_C2C_precorrection_old_();
                    
                    % assemble Coil2Coil precorrection
                    obj.assemble_C2C_precorrection_();
                    
                    % update prev coil to current one
                    obj.setgetPrevCoil(obj.coil);
                    obj.setget_Prev_Zc_prec(obj.Z_C2C);
                end
                
                % form the resulting precorrection
                obj.Z_tot = [obj.Z_C2C, obj.Z_C2B.Nop.';
                    -obj.Z_C2B.Nop, sparse(N_b,N_b)];
            end
                        

                       
           
        end
        
    end
    
% ====================================================================== %

    methods (Access = private)
        
        function assemble_C2B_precorrection_(obj)
            % private method assemble_C2B_precorrection_() precorrects the
            % near interactions between triangle basis functions 
            
            % compute direct C2B interactions
            [Zc2b_Nop_direct, Zc2b_Kop_direct] = obj.assemble_direct_B2C_();
            
            % compute voxelized C2B interactions
            [Zc2b_Nop_voxel, Zc2b_Kop_voxel]  = obj.assemble_voxel_B2C_();
            
            % update the Z_c2b precorrection
            obj.Z_C2B.Nop = sparse(Zc2b_Nop_direct - Zc2b_Nop_voxel);  
            obj.Z_C2B.Kop = sparse(Zc2b_Kop_direct - Zc2b_Kop_voxel);
            
        end
        
        % --------------------------------------------------------------- %
        
        function assemble_C2C_precorrection_(obj)
            % private method assemble_C2C_precorrection_() precorrects the
            % near interactions between voxelized representations of
            % triangles with the directly evaluated ones 
            
            tic;
            
            % assemble Coil2Coil interactions via voxelized representations
            Z_voxel_C2C_near = obj.assemble_voxel_C2C_near_();
            
            t_c2c_vox = toc;
            
            fprintf("Time spent to compute C2C vox precorrection is: %f \n", t_c2c_vox);
            
            
            % assemble direct Coil2Coil interactions 
            Z_direct_C2C_near = obj.assemble_direct_C2C_near_();
            
            % compute Coil2Coil precorrection
            obj.Z_C2C = sparse(Z_direct_C2C_near - Z_voxel_C2C_near);
            
        end
        
        % --------------------------------------------------------------- %
        
        function [Zc2b_Nop_direct, Zc2b_Kop_direct] = assemble_direct_B2C_(obj)
            % method assemlbe_direct_C2B_() computes direct Coil 2 Body
            % interactions 
            
            % compute Coil2Body interactions directly
            [Zc2b_Nop_direct, Zc2b_Kop_direct] = src_precorrection.assemble_direct_B2C(obj.coil, obj.scatterer,...
                                                                 obj.projector, obj.task_settings,...
                                                                 obj.dims, obj.freq);
                                                             
            % scale the Coil2Body interactions with a voxel volume
            Zc2b_Nop_direct = obj.scatterer.dom_ext.res.^3 * Zc2b_Nop_direct;
            Zc2b_Kop_direct = obj.scatterer.dom_ext.res.^3 * Zc2b_Kop_direct;

            
        end
        
        % --------------------------------------------------------------- %
        
        function [Zc2b_Nop, Zc2b_Kop] = assemble_voxel_B2C_(obj)
            % method assemble_voxel_C2B_() computes Coil2Body interactions 
            % using voxelized representation
            
            % get EM constants
            emu      = src_utils.EM_utils(obj.freq);
            
            % verify if gpu flag is set
            gpu_flag = obj.task_settings.general.GPU_flag;
                        
            % assemble B2C interactions with voxelized representation 
            [Zc2b_Nop, Zc2b_Kop]  = src_precorrection.assemble_voxel_B2C(obj.scatterer, obj.projector,...
                                                               obj.mvp, obj.dims, gpu_flag);
            
            % scale resulting matrix                                               
            Zc2b_Nop = 1 / emu.ce * Zc2b_Nop;  
        end
        
        % --------------------------------------------------------------- %
        
        function Z_voxel_C2C_full = assemble_voxel_C2C_near_(obj)
            % method assemble_voxel_C2C_() computes Coil2Coil interactions
            % using voxelized representation
            
            % verify if gpu flag is set
            gpu_flag = obj.task_settings.general.GPU_flag;
                
            Z_voxel_C2C_full = src_precorrection.assemble_voxel_C2C_near(obj.mvp, obj.projector,...
                                                                         obj.dims, obj.freq, gpu_flag);
            
            
        end
        
        % --------------------------------------------------------------- %
        
        
        function assemble_C2C_precorrection_old_(obj)
            % private method assemble_C2C_precorrection_() precorrects the
            % near interactions between voxelized representations of
            % triangles with the directly evaluated ones 
            
            % verify if gpu flag is set
            gpu_flag = obj.task_settings.general.GPU_flag;
            
            % assemble Coil2Coil interactions via voxelized representations
            Z_voxel_C2C_full = src_precorrection.assemble_voxel_C2C(obj.mvp, obj.dims,...
                                                                    obj.freq, gpu_flag);
            
            % compute Coil2Coil precorrection
            obj.Z_C2C = obj.Zc - Z_voxel_C2C_full;
            
        end
        
        % --------------------------------------------------------------- %
        
        function Zc_near = assemble_direct_C2C_near_(obj)
            % private method assemble_direct_C2C_(obj) masks out distant
            % interactions (to be used at precorrection step
            
            % allocate memory for near direct C2C interactions
            Zc_near = zeros(size(obj.Zc));
            
            % get a list of near interactions
            c2c_near_list = obj.projector.C2C_near_list;
            
            % get near direct interactions
            for i = 1:obj.dims.N_sie
                Zc_near(i, c2c_near_list(num2str(i))) = obj.Zc(i, c2c_near_list(num2str(i)));
            end

        end
    end
        % --------------------------------------------------------------- %
    methods(Static)
        
        function prev_Zc_prec = setget_Prev_Zc_prec(Zc_precorrection)
            
            persistent Prev_Zc_prec;
            
            if nargin
                Prev_Zc_prec = Zc_precorrection;
            end
            
            prev_Zc_prec = Prev_Zc_prec;
        end
        % --------------------------------------------------------------- %

        function prevCoil = setgetPrevCoil(coil)
            
            % define persistent variable
            persistent PrevCoil;
            
            % update PrevCoil
            if nargin
                PrevCoil = coil;
            end
            
            % return previous coil
            prevCoil = PrevCoil;
            
        end
        
    end
    
end