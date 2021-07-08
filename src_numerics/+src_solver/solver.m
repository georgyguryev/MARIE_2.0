classdef solver < handle 
        
    properties (SetAccess = immutable)
        
        precond
        
        precond_type
                
        tolerance
        
        max_iter
                
        mvp
        
        dims
        
    end
    
    % ================================================================== %
    
    properties (Access = protected)
        % 
        iterative_solver
    end
    
    % ================================================================== %

    
    methods
        
        function  obj = solver(mvp, precond, max_iter, tolerance, precond_type, dims)
            % constructor of class solver
            
            % get dims object (is empty for decoupled solver)
            obj.dims = dims;
            
            % set up solver parameters
            obj.mvp     = mvp;
            obj.precond = precond;
            
            % set preconditioner and solver types 
            obj.max_iter     = max_iter;
            obj.tolerance    = tolerance;
            obj.precond_type = precond_type;
        end
        
        
        
        function [Jout, rel_res, res_vec] = run(obj, rhs, ini_guess)
            % public method run 
            
            % run iterative solver
            [Jout,~,rel_res,~,res_vec] = obj.iterative_solver(rhs, ini_guess);

        end
        
    end
        
    
    % ================================================================== %
    
    methods (Access = protected, Abstract)
        
        setup_solver_(obj);
         
    end
    
end