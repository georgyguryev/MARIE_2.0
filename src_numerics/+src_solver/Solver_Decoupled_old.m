classdef Solver_Decoupled_old < src_solver.Solver_Base
    
% ====================================================================== %
    
    properties
        outer_iterative_solver
    end
    
% ====================================================================== %
    
    methods
        
        function obj = Solver_Decoupled(mvp, preconditioner, task_settings)
            
            % call constructor of the base class
            obj = obj@src_solver.Solver_Base(mvp, preconditioner, task_settings);
        end
        
        % ------------------------------------------------------------- %
        
        function setup_solver(obj)
            % method setup_solver()
            
            % get solver 
            vie_solver_type = obj.task_settings.vie.Iterative_method;
            vsie_solver_type = obj.task_settings.vsie.Iterative_method;
            
            precond_type = obj.task_settings.vsie.Precond_type;
            
            max_iter_inner = obj.task_settings.vie.MaxIter;
            max_iter_outer = obj.task_settings.vsie.MaxIter;
           
            % get handle to the iterative solvers
            obj.inner_iterative_solver = obj.setup_iterative_solvers_(vie_solver_type, max_iter_inner, precond_type, obj.task_settings);
            obj.outer_iterative_solver = obj.setup_iterative_solvers_(vsie_solver_type, max_iter_outer, precond_type, obj.task_settings);
        end
        
        % ------------------------------------------------------------- %
        
        function [Jcb, rel_res, res_vec] = run(obj,rhs, ini_guess)
            
            % tolerance 
            tol = obj.task_settings.vsie.Tolerance;

            % set up inner VIE solver 
            decoupled_inner_solver = obj.inner_iterative_solver(@(Jin)obj.mvp.VIE(Jin), obj.preconditioner.vie(), tol, []);
            
            % update final VSIE mvp with inner solver
            obj.mvp.update_final_mvp_decoupled_implicit(decoupled_inner_solver);
            
            % set up outer VSIE solver
            decoupled_outer_solver = obj.outer_iterative_solver(@(Jin)obj.mvp.VSIE(Jin),...
                                                                obj.preconditioner.sie(), tol, []);

            % find solution for surface current densities Jc
            [Jc, rel_res, res_vec] = decoupled_outer_solver.run(rhs, ini_guess);
                              
            % compute fields, produced by the surface currents in the
            % scatterer
            Vbc = obj.mvp.c2b_implicit_coupling(Jc);
            
            % find resulting polarization currents
            [Jb,~,~] = decoupled_inner_solver.run(Vbc,[]);
            
            % concatenate solution
            Jcb = [Jc; Jb];
            
        end
        
        % ------------------------------------------------------------- %

    end
    
% ====================================================================== %

end