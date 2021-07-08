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
        
        Zc_inv
                
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
        
        function load_basis(obj)
            
            path = obj.task_settings.basis.Path;
            file = obj.task_settings.basis.Filename;
            
            basis_fname = fullfile(path,file);
            
            load(basis_fname, 'basis');
            
            obj.operator.U_N = basis.U_N;
            obj.operator.U_K = basis.U_K;
            obj.operator.X_N = basis.X_N;
            obj.operator.W_N = basis.W_N;
            obj.operator.alpha_N = basis.alpha_N;
            obj.operator.alpha_K = basis.alpha_K;
            obj.domain = basis.domain;
            obj.index  = basis.index;
        end
        
        
        
        function save_basis(obj)
            
            path = obj.task_settings.basis.Path;
            file = obj.task_settings.basis.Filename;
            
            basis_fname = fullfile(path,file);
            
            index_vie.S_1d = obj.scatterer.index_vie.S_1d;
            obj.domain = obj.scatterer.dom_vie;
            obj.index  = index_vie;
            
            basis = struct('X_N', obj.operator.X_N, 'W_N', obj.operator.W_N,...
                'alpha_N', obj.operator.alpha_N, 'alpha_K', obj.operator.alpha_K,...
                'U_N', obj.operator.U_N, 'U_K', obj.operator.U_K,...
                'domain', obj.domain, 'index', obj.index);
            
            save(basis_fname, 'basis', '-v7.3');
            
            clear basis;
            
        end
        
        
        function assemble_basis_coupling(obj)
            
            N_scat = obj.dims.N_scat;
            
            path = obj.task_settings.basis.Path;
            file = obj.task_settings.basis.Filename;
            
            basis_fname = fullfile(path,file);
           
            tic;
            % check if basis exists 
            if 2 ~= exist(basis_fname,'file')
                obj.construct_basis();
            else
                obj.load_basis();
            end
                       
            [idx, ~] = src_scatterer.get_scat2basis_indices(obj.domain, obj.index,...
                                                            obj.scatterer, obj.freq);
            toc;
            
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
            
            b = obj.operator.M_q2ql * (obj.operator.U_N * (obj.operator.alpha_N * obj.operator.Jc_ini)); 
            
            tic;
            
            
            [Q, S_U, V_hat] = truncated_svd(obj.operator.U_N, 1e-6);
            toc;
            tic;
            W = S_U * V_hat' * obj.operator.X_N * obj.operator.alpha_N.' * conj(V_hat) * S_U;

            [Q_w,S_w, V_w] = truncated_svd(W, 1e-20);
            
            U = Q * Q_w;
            alpha = U' * b;
            c = sqrt(S_w * alpha);
            [c_sort, J] = sort(c, 'descend');
            
            basis_rank = find_max_rank(c_sort, obj.task_settings.vsie.Tolerance);
            
            J_unique = unique(reshape(J(1:basis_rank,:), [],1));
            
            UU_N = Q * Q_w * S_w;
            V_N  = V_w' * Q.';
            
            obj.operator.UU_N = UU_N(:,J_unique);
            obj.operator.V_N = V_N(J_unique,:);
            
            toc;
            
            %% check if operators shouls be moved to GPU
            if obj.task_settings.basis.GPU_flag
                obj.operator.U_N     = src_utils.to_GPU(obj.operator.U_N, 0);
                
                obj.operator.UU_N    = src_utils.to_GPU(obj.operator.UU_N, 0);
                obj.operator.V_N     = src_utils.to_GPU(obj.operator.V_N, 0);
%                 obj.operator.U_K     = src_utils.to_GPU(obj.operator.U_K, 0);
                obj.operator.alpha_N = src_utils.to_GPU(obj.operator.alpha_N, 0);
%                 obj.operator.alpha_K = src_utils.to_GPU(obj.operator.aplha_K, 0);
                obj.operator.X_N     = src_utils.to_GPU(obj.operator.X_N, 0);
                obj.operator.W_N     = src_utils.to_GPU(obj.operator.W_N, 0);
                obj.operator.M_q2ql  = src_utils.to_GPU(obj.operator.M_q2ql, 0);
            end
          
        end        
    end
    
    methods (Abstract)
        
        % declare method for basis assembly
        construct_basis(obj)
                            
    end
    
end