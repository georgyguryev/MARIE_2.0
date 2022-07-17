classdef Solver_Tissue_Implicit < src_solver.Solver_Base
    
% ====================================================================== %
    
    properties
        outer_iterative_solver
    end
    
% ====================================================================== %
    
    methods
        
        function obj = Solver_Tissue_Implicit(mvp, preconditioner, task_settings)
            
            % call constructor of the base class
            obj = obj@src_solver.Solver_Base(mvp, preconditioner, task_settings);
        end
        
        % ------------------------------------------------------------- %
        
        function setup_solver(obj)
            % method setup_solver()
            
%             % get solver 
%             solver_type = obj.task_settings.vsie.Iterative_method;
%             
%             max_iter     = obj.task_settings.vsie.MaxIter;
%             precond_type = obj.task_settings.vsie.Precond_type;
%            
%             % get handle to the iterative solvers
%             obj.inner_iterative_solver = obj.setup_iterative_solvers_(solver_type, max_iter, precond_type, obj.task_settings);
            
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
        
        function [Jb, rel_res, res_vec] = run(obj,rhs, ini_guess, feed_port)
            
            % tolerance 
            dims    = obj.mvp.dims;
            cur_rhs = rhs(:,feed_port);   

%             % get tolerance
%             tolerance  = obj.task_settings.vsie.Tolerance;

            % get preconditioner
%             preconditioner = obj.preconditioner.vie();
%             preconditioner = obj.preconditioner.vie();
%             
%             decoupled_solver = obj.inner_iterative_solver(@(Jin)obj.mvp.VSIE(Jin), preconditioner, tolerance, dims);
% 
%             % run iterative solver
%             [Jc, rel_res, res_vec] = decoupled_solver.run(cur_rhs, ini_guess);           
%             
%             % compute surface currents 
%             Jb  = obj.mvp.VIE(Jc, feed_port);
            
            % concatenate solution
%             Jcb = [Jc; Jb];


            % tolerance 
            tol = obj.task_settings.vsie.Tolerance;

            % set up inner VIE solver 
            decoupled_inner_solver = obj.inner_iterative_solver(@(Jin)obj.mvp.VIE(Jin), obj.preconditioner.vie(), tol, dims);
            
            % update final VSIE mvp with inner solver
            obj.mvp.update_final_mvp_decoupled_implicit(decoupled_inner_solver);
            
            % set up outer VSIE solver
            decoupled_outer_solver = obj.outer_iterative_solver(@(Jin)obj.mvp.VSIE(Jin), obj.preconditioner.sie(), tol, dims);

            % find solution for surface current densities Jc
            [Jc, rel_res, res_vec] = decoupled_outer_solver.run(cur_rhs, ini_guess);
                              
            % compute fields, produced by the surface currents in the
            % scatterer
            Vbc = obj.mvp.c2b_implicit_coupling(Jc);
            
            % find resulting polarization currents
            [Jb,~,~] = decoupled_inner_solver.run(Vbc,[]);
            
            % concatenate solution
%             Jcb = [Jc; Jb];
            
        end
        
        % ------------------------------------------------------------- %

    end
    
% ====================================================================== %

end