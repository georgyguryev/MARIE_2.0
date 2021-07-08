classdef Cross_approx < Sampling_approx_Base
    
    
    properties
        
        fV
        fA_hat
        
    end
    
    
    methods
        
        function obj = Cross_approx(task_settings, scatterer, coil, dims, freq)
        % constructor of a cross_approximation object
            
            obj = obj@Sampling_approx_Base(task_settings, scatterer, coil, dims, freq);
          
%             obj.fV = @(rows, op_type) get_V_bb(rows, obj.task_settings,...
%                                                     obj.scatterer, obj.coil,...
%                                                     obj.dims, freq, op_type);
                                                
            obj.fV = @(rows, op_type) get_V(rows, obj.task_settings,...
                                                    obj.scatterer, obj.coil,...
                                                    obj.dims, freq, op_type); 
                                                
%             obj.fA_hat = @(rows, cols, op_type) get_A_hat_bb(rows, cols, obj.task_settings,...
%                                                                 obj.scatterer, obj.coil,...
%                                                                 obj.dims, freq, op_type); 

           obj.fA_hat = @(rows, cols, op_type) get_A_hat(rows, cols, obj.task_settings,...
                                                                obj.scatterer, obj.coil,...
                                                                obj.dims, freq, op_type); 
        end
        
        
        
        function [U,V] = compute_approximation(obj, op_type)
            % initialize row and column sampling methods
            
            % get dimensions 
            N = obj.dims.N_sie;            
            M = obj.dims.ql * obj.dims.N_scat;
            
            % cross-related settings
            tol      = obj.task_settings.basis.Tolerance;
            Max_iter = obj.task_settings.basis.Max_iter_cross;           
            
            % specify operator type ('N' or 'K';  'NK' otherwise)
            fU_op = @(cols) obj.fU(cols, op_type);
            fV_op = @(rows) obj.fV(rows, op_type);
            fA_hat_op = @(rows, cols) obj.fA_hat(rows, cols, op_type);
            
            profile on;                                                      

            idxS = obj.scatterer.index_vie.S_1d;
            
            % run cross 2d
            [U,V] = cross_cheb_2d(fU_op, fV_op, fA_hat_op, obj.fTerm,...
                                    idxS, N, obj.dims, Max_iter, tol);
                                
            profile off;
            profile viewer;
        end      
        
    end
    
    
end
