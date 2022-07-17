classdef Coupling_Basis_Cross < Coupling_Basis_Base  
    
    properties
        
        cross
        
        % store left singular vectors for N and K operators
        V_N
        V_K
        
    end
    
    
    methods
        
        function obj = Coupling_Basis_Cross(task_settings, scatterer, coil, operator, dims, freq, Zc_inv, Z_L, F, rhs_cp)
            
            obj = obj@Coupling_Basis_Base(task_settings, scatterer, coil, operator, dims, freq);
            
            % set up coil-specific matrices if basis has to be assembled
            if nargin > 6
                obj.Zc_inv = Zc_inv;
                obj.Z_L    = Z_L;
                obj.F      = F;
                obj.rhs_cp = rhs_cp;
                Jini       = Zc_inv * F;
                
                obj.Zc_inv_hat = Zc_inv + Jini / (inv(Z_L) - F.' * Jini) * Jini.';
            end
            
        end
        
        % --------------------------------------------------------------- %
        
        function  construct_perturbation_basis(obj)
            
            % instantiate cross approximation 
            obj.cross = Cross_approx(obj.task_settings, obj.scatterer, obj.coil, obj.dims, obj.freq);
            
            profile on;
            % compute N and K crosses
            [U_N, alpha_N] = obj.cross.compute_approximation('N');
            profile off;
            profile viewer;
            keyboard; 

%             profile on;            
            [U_K, alpha_K] = obj.cross.compute_approximation('K');
%             profile off;
%             profile viewer;
            
            X_N = alpha_N * obj.Zc_inv;
            X_K = alpha_K * obj.Zc_inv;

            
            W_N = X_N * alpha_N.';
            W_K = X_K * alpha_N.';
            
            N_scat    = obj.dims.N_scat;            
            
            b_N = U_N * (alpha_N * (obj.Zc_inv * obj.F)); 
            b_K = U_K * (alpha_K * (obj.Zc_inv * obj.F));
            b_Icu = -obj.F.' * obj.Zc_inv * obj.F;
            
            
            % compress basis from rank(Zbc) to rank(Zbb^p)
            [Qw, Sw, Vw]    = truncated_svd(W_N, 1e-20);
            [Qwk, Swk, Vwk] = truncated_svd(W_K, obj.task_settings.basis.Tolerance);
            
            Q_N = U_N * Qw; 
            Q_K = U_K * Qwk;
            
            % determine basis expansion coefficients with the largest
            % weight
%             alpha = Q_N' * b_N;
            
%             c = sqrt(Sw * alpha);
%             [c_sort, J] = sort(c, 'descend');
            
%             basis_rank = find_max_rank(c_sort, obj.task_settings.basis.Tolerance);

%             J_unique = unique(reshape(J(1:basis_rank,:), [],1));

            basis_rank = find_max_rank(diag(Sw), obj.task_settings.basis.Tolerance);
            
            V_N = Sw * Vw' * U_N.';
            V_K = Swk * Vwk' * U_N.';
            
%             X_cs = obj.F.' * obj.Zc_inv_hat * alpha_N.' * U_N.';
            X_cu = obj.F.' * obj.Zc_inv * alpha_N.' * U_N.';
            
            
            obj.operator.U_N = Q_N(:,1:basis_rank);
            obj.operator.U_K = Q_K;
            obj.operator.V_N = V_N(1:basis_rank,:);
            obj.operator.V_K = V_K;
            obj.operator.X_cu = X_cu;
            obj.operator.b_N = b_N;
            obj.operator.b_K = b_K;
            obj.operator.b_Icu = b_Icu; 
            obj.operator.Jc_ini = obj.Zc_inv * obj.F;
            
            % save basis
            obj.save_basis();
        end
        
        % --------------------------------------------------------------------------------%
        
        function  construct_coupled_basis(obj)
            
            % instantiate cross approximation 
            obj.cross = Cross_approx(obj.task_settings, obj.scatterer, obj.coil, obj.dims, obj.freq);
            profile on;
            % compute N and K crosses
            [U_N, alpha_N] = obj.cross.compute_approximation('N');
            [U_K, alpha_K] = obj.cross.compute_approximation('K');
            
            keyboard;
                       
            b_N = U_N * (alpha_N * obj.operator.Jc_ini); 
            b_K = U_K * (alpha_K * (obj.Zc_inv * obj.F));
            
            b_Icu = -obj.F.' * obj.Zc_inv * obj.F;
            X_cu = obj.F.' * obj.Zc_inv * alpha_N.' * U_N.';
            
            X_K = alpha_K * obj.Zc_inv;
            
            W_K = X_K * alpha_N.';
            
            [Qwk, Swk, Vwk] = truncated_svd(W_K, obj.task_settings.basis.Tolerance);
            
            Q_K = U_K * Qwk;
            V_K = Swk * Vwk' * U_N.';

            
            obj.operator.U_N = U_N;
            obj.operator.U_K = Q_K;
            obj.operator.V_N = alpha_N;
            obj.operator.V_K = V_K;
            obj.operator.X_cu = X_cu;
            obj.operator.b_N = b_N;
            obj.operator.b_K = b_K;
            obj.operator.b_Icu = b_Icu; 
            obj.operator.Jc_ini = obj.Zc_inv * obj.F;
            
            % save basis
            obj.save_basis();
        end
                
    end
   
end