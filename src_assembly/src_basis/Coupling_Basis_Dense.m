classdef Coupling_Basis_Dense < Coupling_Basis_Base  
    
    properties
        
        I_UN
        I_UK
        fV
    end
    
    
    methods
        
        function obj = Coupling_Basis_Dense(task_settings, scatterer, coil, operator, dims, Zc_inv, freq)
            
            
            obj = obj@Coupling_Basis_Base(task_settings, scatterer, coil, operator, dims, freq);
            
            obj.fV = @(rows, op_type) get_V(rows, obj.task_settings,...
                                            obj.scatterer, obj.coil,...
                                            obj.dims, freq, op_type);                            
            obj.Zc_inv = Zc_inv;
                            
        end
        
        % --------------------------------------------------------------- %
        
        function  construct_basis(obj)
            
            % get singular vectors
            obj.construct_basis_();
            
            % save basis
            obj.save_basis();
        end
     
    end
    
 % ---------------------------------------------------------------------- %   
    methods (Access = private)
                
        function construct_basis_(obj)
            
            tol = obj.task_settings.basis.Tolerance;

            % instantiate cross approximation 
            [Z_N, Z_K] = src_coupling.assemble_coupling_matrices(obj.task_settings, obj.scatterer,...
                                                                 obj.coil,  obj.dims, obj.freq,[],[],0);
            % compute singular vectors
            [U_N,S_N,V_N] = svd(Z_N, 'econ');
            [U_K,S_K,V_K] = svd(Z_K, 'econ');
            
            clear Z_N Z_K;
            
            %  find rank and truncate singular vectors
            rank_N = find_rank_from_S(diag(S_N), tol);
            rank_K = find_rank_from_S(diag(S_K), tol);
            
            obj.operator.U_N = U_N(:,1:rank_N);
            obj.operator.U_K = U_K(:,1:rank_K);
            obj.operator.alpha_N = S_N(1:rank_N, 1:rank_N) * V_N(:,1:rank_N)';
            obj.operator.alpha_K = S_K(1:rank_K, 1:rank_K) * V_K(:,1:rank_K)';
            obj.operator.X_N = obj.operator.alpha_N * obj.Zc_inv;
            obj.operator.W_N = obj.operator.X_N * obj.operator.alpha_N.';
            
            clear U_N U_K alpha_N alpha_K X_N W_N;
            
        end
        
        % --------------------------------------------------------------- %
       
        function construct_DEIM_mask_(obj)
            
            % check if DEIM++ is available
            if (3 == exist('mexDeim','file'))
                [~, obj.Pin_N] = mexDeim(obj.U_N);
                [~, obj.Pin_K] = mexDeim(obj.U_K);
            else
                [~, obj.Pin_N] = src_numeric.deim(obj.U_N);
                [~, obj.Pin_K] = src_numeric.deim(obj.U_K);
            end
            
        end
        
        % --------------------------------------------------------------- %
       
        function postprocess_basis_(obj)
            %% postprocessing
            
            obj.I_UN = maxvol(obj.operator.U_N, 1e-7, min(size(obj.operator.U_N)));
            obj.I_UK = maxvol(obj.operator.U_K, 1e-7, min(size(obj.operator.U_K)));

            UN_hat   = obj.operator.U_N(obj.I_UN,:);
            UK_hat   = obj.operator.U_K(obj.I_UK,:); 
            
            ZbcN_hat = obj.fV(obj.I_UN, 'N');
            ZbcK_hat = obj.fV(obj.I_UK, 'K');
            
            obj.operator.alpha_N = UN_hat \ ZbcN_hat;
            obj.operator.alpha_K = UK_hat \ ZbcK_hat;
            obj.operator.X_N = obj.operator.alpha_N * obj.Zc_inv;
            obj.operator.W_N = obj.operator.X_N * obj.operator.alpha_N.';

        end
        
        
        
    end
   
end