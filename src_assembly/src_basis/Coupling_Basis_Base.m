classdef Coupling_Basis_Base < handle 
    
    properties
        
        % input parameters
        task_settings
        scatterer
        coil
        operator
        dims
        freq
        
        domain
        index
        
        % store left singular vectors for N and K operators
        U_N
        U_K
        X_cu
        X_cs
        
        b_N
        b_K
        b_Ics
        b_Icu
        
        Zc_inv
        Zc_inv_hat
        Z_L  
        
        F
        rhs_cp
    end
    
    
    methods
        
        function obj = Coupling_Basis_Base(task_settings, scatterer, coil, operator, dims, freq)
            
            % constructor of class Coupling_Basis_Base
            obj.task_settings = task_settings; 
            obj.scatterer     = scatterer;
            obj.operator      = operator;
            obj.coil          = coil;
            obj.dims          = dims;
            obj.freq          = freq;
        end
        
        % --------------------------------------------------------------- %

        function load_basis(obj)
            
            path = obj.task_settings.basis.Path;
            file = src_utils.prefix_basis_filename(obj.task_settings);
            
%             % label basis type with prefix
%             switch solver_mode
%                 case 'Coupled'
%                     file = strcat('Coupled_basis_', file);
%                 case 'Decoupled'
%                     file = strcat('Perturbation_basis_', file);
%             end
            
            basis_fname = fullfile(path,file);
            
            load(basis_fname, 'basis');
            
            obj.operator.U_N = basis.U_N;
            obj.operator.U_K = basis.U_K;
            obj.operator.b_N = basis.b_N;
            obj.operator.b_K = basis.b_K;
            obj.operator.b_Ics = basis.b_Ics;
            obj.operator.b_Icu = basis.b_Icu;           
            
            obj.operator.X_cs = basis.X_cs;
            obj.operator.X_cu = basis.X_cu;

            obj.operator.V_N = basis.V_N;
            obj.operator.V_K = basis.V_K;
            
            obj.operator.Jc_ini = basis.Jc_ini; 
            
            obj.domain = basis.domain;
            obj.index  = basis.index;
        end
        
        % --------------------------------------------------------------- %
        
        
        function save_basis(obj)
            
            path = obj.task_settings.basis.Path;
            file = src_utils.prefix_basis_filename(obj.task_settings);
            
            basis_fname = fullfile(path,file);
            
            index_vie.S_1d = obj.scatterer.index_vie.S_1d;
            obj.domain = obj.scatterer.dom_vie;
            obj.index  = index_vie;
            
            basis = struct('X_cs', obj.operator.X_cs, 'X_cu', obj.operator.X_cu,...
                'V_N', obj.operator.V_N, 'V_K', obj.operator.V_K,...
                'U_N', obj.operator.U_N, 'U_K', obj.operator.U_K,...,
                'b_N', obj.operator.b_N, 'b_K', obj.operator.b_K,...
                'b_Ics', obj.operator.b_Ics, 'b_Icu', obj.operator.b_Icu,...
                'domain', obj.domain, 'index', obj.index, 'Jc_ini', obj.operator.Jc_ini);
            
            save(basis_fname, 'basis', '-v7.3');
            
            clear basis;
            
        end
        
        % --------------------------------------------------------------- %
        
        
        
        function assemble_basis_coupling(obj)

            path = obj.task_settings.basis.Path;
            file = src_utils.prefix_basis_filename(obj.task_settings);
            solver_mode = lower(obj.task_settings.vsie.Solver_mode);
            
            basis_fname = fullfile(path,file);
           
            tic;
            % check if basis exists 
            if 2 ~= exist(basis_fname,'file')
                profile on;
                
                switch solver_mode
                    case 'explicit'
                        obj.construct_coupled_basis();
                    case 'coil_implicit'
                        obj.construct_perturbation_basis();
                    case 'tissue_implicit'
                        obj.construct_coupled_basis();
                end
                
                profile off;
                profile viewer;
                keyboard;
            else
                obj.load_basis();
            end
                       
            [idx, ~] = src_scatterer.get_scat2basis_indices(obj.domain, obj.index,...
                                                            obj.scatterer, obj.freq);
            toc;
            N_scat = obj.dims.N_scat;
            Ndofs_pwx = obj.dims.ql * N_scat;
              
            if size(obj.operator.U_N,1) == Ndofs_pwx 
                obj.operator.M_q2ql = speye(Ndofs_pwx);
            else
                idx_q = find(kron(ones(obj.dims.q,1),idx));
                obj.operator.U_N = obj.operator.U_N(idx_q, :);
                obj.operator.U_K = obj.operator.U_K(idx_q, :);
                
                idx_Sl                     = sparse(N_scat * obj.dims.l, N_scat);
                idx_Sl(1:N_scat, 1:N_scat) = speye(N_scat);
                obj.operator.M_q2ql = kron(speye(obj.dims.q), idx_Sl);
            end
            
            %% check if operators shouls be moved to GPU
            if obj.task_settings.basis.GPU_flag
                obj.operator.U_N    = src_utils.to_GPU(obj.operator.U_N, 0);
                obj.operator.V_N    = src_utils.to_GPU(obj.operator.V_N, 0);
            end
          
        end     
        
        % --------------------------------------------------------------- %
        
       
        
    end
    
    methods (Abstract)
        
        % declare method for basis assembly
        construct_perturbation_basis(obj)
        
        construct_coupled_basis(obj)
                            
    end
    
end