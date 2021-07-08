classdef Solver_Base < handle

% ====================================================================== %

    properties %(Access = protected)
        mvp                     % matrix-vector product
        preconditioner          % preconditioner
        task_settings           % task settings
        inner_iterative_solver  % iterative solver
    end
    
% ====================================================================== %
    
    methods
        
        function obj = Solver_Base(mvp, preconditioner, task_settings)
            
            obj.mvp            = mvp;
            obj.preconditioner = preconditioner;
            obj.task_settings  = task_settings;
            
        end
    end
    
% ====================================================================== %
    
    methods (Access = protected)
        
        function iterative_solver = setup_iterative_solvers_(~, iter_solver_type,  max_iter, precond_type,task_settings)
            
            switch iter_solver_type
                case 'GMRES_DR'
                    iterative_solver = @(mvp, precond, tolerance, dims) src_solver.gmres_DR_solver(@(Jin) mvp(Jin), precond, max_iter,...
                                         tolerance, precond_type, task_settings, dims);
                                     
                case 'GMRES'
                    iterative_solver = @(mvp, precond,tolerance, dims) src_solver.gmres_solver(@(Jin) mvp(Jin), precond,  max_iter,...
                                         tolerance, precond_type, task_settings, dims);
                                     
                case 'BICGSTAB'
                    iterative_solver = @(mvp, precond, tolerance, dims) src_solver.bicgstab_solver(@(Jin) mvp(Jin), precond,  max_iter,...
                                        tolerance, precond_type, dims);
                                    
                case 'TFQMR'
                    iterative_solver = @(mvp, precond, tolerance, dims)src_solver.tfqmr_solver(@(Jin) mvp(Jin), precond,  max_iter,...
                                        tolerance, precond_type, dims);
                                    
            end
        end
    end
    
end