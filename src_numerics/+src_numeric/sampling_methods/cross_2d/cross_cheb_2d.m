function  [X, W] = cross_cheb_2d(fU, fV, fA_hat, f_term_criteria, idxS, N, dims,Max_iter, tol)
    
    n0 = min(dims.vie(1:3));
    
    maxvol_tol = tol;
    max_rank = 0;
    rank = 100; %round(dims.N_sie / 10);
    
    I = sample_rows([], dims, rank, n0, 1, idxS);
    J = sample_cols([], [], rank, N, maxvol_tol);
    
    AJ = fU(J);
    AI = fV(I);

    % orthogonalize  
    [U,~] = qr(AJ,0);
    [V,~] = qr(AI',0); 
    
    S_new = 0;
    old_rank = 0;
    
    for iter = 2:Max_iter
        
        [J_hat, rank_hat] = sample_cols(V, J, rank, [], maxvol_tol);
        I_hat             = sample_rows(I, dims, rank_hat, n0, iter, idxS);
        
        if rank_hat == 0
            break;
        end
        
        I = union(I, I_hat);
        J = union(J, J_hat);
        
        % sample selected columns and rows       
        U_aux = fU(J_hat);
        V_aux = fV(I_hat);
        
        % form a new submatrix
        U = [U, U_aux];
        V = [V, V_aux'];
                
        [U,~] = qr(U,0);
        [V,~] = qr(V,0);
        
        U_hat = U(I,:);
        V_hat = V(J,:);
        A_hat = fA_hat(I, J);
        
        [Qu_hat,Su_hat,Vu_hat] = svd(U_hat, 'econ'); 
        [Qv_hat,Sv_hat,Vv_hat] = svd(V_hat', 'econ');
        
        M = Su_hat \ (Qu_hat' * A_hat * Vv_hat) / Sv_hat;
        
        [Qm,Sm,Vm] = truncated_svd(M, tol);        
        
        Qu = Vu_hat * Qm * Sm;
        Qv = Vm' * Qv_hat';
        
        S_old = S_new; 
        S_new = diag(Sm);
        rank = find_rank_from_S(S_new, tol);
        
        max_rank = max(max_rank,rank);
        
        X = U * Qu;
        W = Qv * V';  
        
        U = U(:,1:max_rank);
        V = V(:,1:max_rank);
              
        % compute norm of difference in singular values (old - new)
        if f_term_criteria(S_old, S_new, old_rank, rank, tol)
            return;
        end
                
       old_rank = rank;
       
    end
end