classdef Solver_Coil_Implicit < src_solver.Solver_Base
    
% ====================================================================== %
    
    properties
        outer_iterative_solver
    end
    
% ====================================================================== %
    
    methods
        
        function obj = Solver_Coil_Implicit(mvp, preconditioner, task_settings)
            
            % call constructor of the base class
            obj = obj@src_solver.Solver_Base(mvp, preconditioner, task_settings);
        end
        
        % ------------------------------------------------------------- %
        
        function setup_solver(obj)
            % method setup_solver()
            
            % get solver 
            solver_type = obj.task_settings.vsie.Iterative_method;
            
            max_iter     = obj.task_settings.vsie.MaxIter;
            precond_type = obj.task_settings.vsie.Precond_type;
           
            % get handle to the iterative solvers
            obj.inner_iterative_solver = obj.setup_iterative_solvers_(solver_type, max_iter, precond_type, obj.task_settings);
        end
        
        % ------------------------------------------------------------- %
        
        function [Jb, rel_res, res_vec] = run(obj,rhs, ini_guess, feed_port)
            
            % tolerance 
            dims    = obj.mvp.dims;
            cur_rhs = rhs(:,feed_port);   

            % get tolerance
            tolerance  = obj.task_settings.vsie.Tolerance;

            % get preconditioner
            preconditioner = obj.preconditioner.vie();
            
            decoupled_solver = obj.inner_iterative_solver(@(Jin)obj.mvp.VSIE(Jin), preconditioner, tolerance, dims);

            % run iterative solver
            [Jb, rel_res, res_vec] = decoupled_solver.run(cur_rhs, ini_guess);           
            
        end
        
        % ------------------------------------------------------------- %

    end
    
% ====================================================================== %

end