                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       classdef Sampling_approx_Base < handle 
  
    properties
        
        task_settings
        scatterer
        coil
        dims
        
        fU 
        fTerm
    end
    
    
    methods
        
        function obj = Sampling_approx_Base(task_settings, scatterer, coil, dims, freq)
            % constructor of a Sampling_approx_Base object
            
            obj.task_settings = task_settings;
            obj.scatterer = scatterer;
            obj.coil = coil;
            obj.dims = dims;
            
%             obj.fU = @(cols, op_type) get_U_bb(cols, obj.task_settings,...
%                                                     obj.scatterer, obj.coil,...
%                                                     obj.dims, freq, op_type);

            obj.fU = @(cols, op_type) get_U(cols, obj.task_settings,...
                                                    obj.scatterer, obj.coil,...
                                                    obj.dims, freq, op_type);                                    
                                                
                                                                        
            obj.fTerm = @(S_old, S, rank_old, rank, tol) obj.Frobenius_termination_criteria(S_old, S, rank_old, rank, tol);
                                                
        end
      
        
        function converged = Frobenius_termination_criteria(~, S_old, S, rank_old, rank, tol)
                       
            rel_accuracy = tol * norm(S(1:rank));
            
            if rank_old < rank
                ss_error     = norm(S_old(1:rank_old) - S(1:rank_old));
%                 ss_error = 0;
                trunc_error  = norm(S(rank_old+1:rank));
            else
                ss_error     = norm(S_old(1:rank) - S(1:rank));
                trunc_error  = 0;
            end
            approx_error = sqrt(ss_error.^2 + trunc_error.^2);
            
            % check termination criteria
             converged = (rel_accuracy > approx_error);
%             converged = (rel_accuracy > trunc_error);
           
        end
        
    end    
    
    
    methods (Abstract)
        
       compute_approximation(obj, op_type)

    end
    
end