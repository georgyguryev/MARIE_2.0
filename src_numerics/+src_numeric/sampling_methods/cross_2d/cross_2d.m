 function  [U, V] = cross_2d(fU, fV, fA_hat, f_term_criteria, M, N, Max_iter, tol)
    
    rank     = 80;
    old_rank = 1;
        
    samp_cols = randi(N, rank, 1);
    samp_rows = randi(M, rank, 1);

    %% ols, op_type, freq)
        
    U = fU(samp_cols);
    V = fV(samp_rows);

    % orthogonalize  
    [U,~] = qr(U,0);
    [V,~] = qr(V.',0);
    V = V.';
    
    S_new = 0;
    
    for iter = 1:Max_iter
        
        I_U = maxvol(U, 0.2, rank);
        I_V = maxvol(V.', 0.2, rank);
        
        %% sample selected columns and rows       
        U_aux = fU(I_V);
        V_aux = fV(I_U);
        
        % form a new submatrix
        U = [U, U_aux];
        V = [V; V_aux];
        
        rank = 2 * rank;
        
        [U,~] = qr(U,0);
        [V,~] = qr(V.',0);
        V = V.';
        
        I_U = maxvol(U, 0.2, rank);
        I_V = maxvol(V.', 0.2, rank);
        
        U_hat = U(I_U,:);
        V_hat = V(:,I_V);
        A_hat = fA_hat(I_U, I_V);
        
        % find G = U_hat^-1 * A_hat * V_hat^-1;
        G = U_hat \ A_hat / V_hat;      
        
        [Qu, S, Qv] = svd(G);
        
        S_old = S_new; 
        S_new = diag(S);
        rank = find_rank_from_S(S_new, tol);
        
        S = diag(S_new(1:rank));
        U = U * Qu(:,1:rank) * S;
        V = Qv(:,1:rank)' * V;
        
        % compute norm of difference in singular values (old - new)
        if f_term_criteria(S_old, S_new, old_rank, rank, tol)
            return;
        end
        
       old_rank = rank;
       
    end
end