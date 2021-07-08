classdef bicgstab_solver < src_solver.solver
    
    % ================================================================== %

    properties (SetAccess = immutable)
     
    end
    
    % ================================================================== %
    
    methods
        
        function obj = bicgstab_solver(mvp, precond,  max_iter, tolerance, precond_type, dims)
            % constructor of class gmres_DR
            
            % call base class constructor 
            obj@src_solver.solver(mvp, precond,  max_iter, tolerance, precond_type, dims);
                        
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
            
            % setup bicgstab iterative method
           obj.iterative_solver = @(rhs, ini_guess)  src_numeric.my_bicgstab(@(Jin) obj.mvp(Jin), rhs, obj.tolerance,...
                        obj.max_iter, Lp, Pp, Up, ini_guess);
         
           
        end
        
    end

end