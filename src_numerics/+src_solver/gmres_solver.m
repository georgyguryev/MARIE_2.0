classdef gmres_solver < src_solver.solver
   
    properties (SetAccess = immutable)
        
        restart
        
    end
    
    % ================================================================== %
    
    methods 
        
        function obj = gmres_solver(mvp, precond, max_iter, tolerance, precond_type,...
                                    task_settings, dims)
                                
            % constructor of class gmres_solver
            
            % call constructor of the base class
            obj@src_solver.solver(mvp, precond,  max_iter, tolerance, precond_type, dims);
            
            % set up a number of iterations before GMRES restart 
            obj.restart = task_settings.vsie.Restart;
            
            % set up solver (select appropriate type of preconditioner) 
            obj.setup_solver_();
        end
        
    end
    
    % ================================================================== %
   
    methods (Access = protected)
        
        function setup_solver_(obj)
            % function setup_solver_() sets up the solver with appropriate
            % preconditioner type
            
            if isempty(obj.precond)
                Lp = []; Up =[]; Pp = [];
            else
                % get preconditioning matrices
                Lp = obj.precond.L;
                Up = obj.precond.U;
                Pp = obj.precond.P;
            end
            
            switch obj.precond_type
                
                case 'Symmetric'
                    obj.iterative_solver = @(rhs, ini_guess) src_numeric.pgmres(@(Jin) obj.mvp(Jin), rhs, obj.restart, obj.tolerance,...
                                                                                obj.max_iter, Lp, Pp, [], Up, [], ini_guess);
                case 'Left'
                    obj.iterative_solver = @(rhs, ini_guess) src_numeric.pgmres(@(Jin) obj.mvp(Jin), rhs, obj.restart, obj.tolerance,...
                                                                                obj.max_iter, Lp, Pp, Up, [], [], ini_guess);
                case 'Right'
                    obj.iterative_solver = @(rhs, ini_guess) src_numeric.pgmres(@(Jin) obj.mvp(Jin), rhs, obj.restart, obj.tolerance,...
                                                                                obj.max_iter, [], [], [], Lp, Up, ini_guess);
            end
            
        end
         
    end
    
end