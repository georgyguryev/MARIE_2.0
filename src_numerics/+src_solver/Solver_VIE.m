classdef Solver_VIE < src_solver.Solver_Base
    
% ====================================================================== %
    
    methods
        
        function obj = Solver_VIE(mvp, preconditioner, task_settings)
            
            % call constructor of the base class
            obj = obj@src_solver.Solver_Base(mvp, preconditioner, task_settings);
        end
        
        % ------------------------------------------------------------- %
        
        function setup_solver(obj)
            % method setup_solver()
            
            % get solver 
            solver_type = obj.task_settings.vie.Iterative_method;
            
            max_iter     = obj.task_settings.vie.MaxIter;
            precond_type = obj.task_settings.vsie.Precond_type;

            % get handle to the iterative solver
            obj.inner_iterative_solver = obj.setup_iterative_solvers_(solver_type, max_iter, precond_type, obj.task_settings);
        end
        
        % ------------------------------------------------------------- %
        
        function [Jcb, rel_res, res_vec] = run(obj,rhs,ini_guess)
            % method run(rhs) 
            
            dims = obj.mvp.dims;
            
            % get tolerance
            tol = obj.task_settings.vie.SolTolerance;
            
            preconditioner = obj.preconditioner.vie();
            
            coupled_solver = obj.inner_iterative_solver(@(Jin)obj.mvp.VIE(Jin), preconditioner, tol, dims);
           
            % run iterative solver
            [Jcb, rel_res, res_vec] = coupled_solver.run(rhs, ini_guess);
                        
        end
        
        % ------------------------------------------------------------- %

    end
    
% ====================================================================== %

end