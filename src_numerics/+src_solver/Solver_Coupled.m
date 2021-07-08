classdef Solver_Coupled < src_solver.Solver_Base
    
% ====================================================================== %
    
    methods
        
        function obj = Solver_Coupled(mvp, preconditioner, task_settings)
            
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

            % get handle to the iterative solver
            obj.inner_iterative_solver = obj.setup_iterative_solvers_(solver_type, max_iter, precond_type, obj.task_settings);
        end
        
        % ------------------------------------------------------------- %
        
        function [Jcb, rel_res, res_vec] = run(obj, rhs, ini_guess,  feed_port)
            % method run(rhs) 
            
            dims    = obj.mvp.dims;
            cur_rhs = rhs(:,feed_port);   
            
            % get tolerance
            tolerance = obj.task_settings.vsie.Tolerance;

            profile off;
            profile viewer;
            

            preconditioner = obj.preconditioner.vsie();
            
            coupled_solver = obj.inner_iterative_solver(@(Jin)obj.mvp.VSIE(Jin), preconditioner, tolerance, dims);
                        
            
            profile on; 
            % run iterative solver
            [Jcb, rel_res, res_vec] = coupled_solver.run(cur_rhs, ini_guess);
            profile off;
            profile viewer;

         
            
            fprintf('Number of iterations: %d; \n', length(res_vec));
            fprintf('Relative residual: %6.4e \n', rel_res);
    
            
        end
        
        % ------------------------------------------------------------- %

    end
    
% ====================================================================== %

end