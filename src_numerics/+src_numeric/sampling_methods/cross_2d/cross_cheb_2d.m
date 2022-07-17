function  [X, W] = cross_cheb_2d(fU, fV, fA_hat, f_term_criteria, idxS, N, dims,Max_iter, tol)
    
%     n_min = min(dims.vie(1:3));
    
    maxvol_tol = 1e-2;
    max_rank = 0;
    rank  = 100; %round(dims.N_sie / 10);
    stall = 0;
    n0 = 2*ceil((rank)^(1/3));   %
%     n0 = min(dims.vie(1:3));

%     I = sample_rows_ql([], dims, rank, n0, 1, idxS);
    I = sample_rows([], dims, rank, n0, 1, idxS);
    J = sample_cols([], [], rank, N, maxvol_tol);
    
    AJ = fU(J);
    AI = fV(I);

    % orthogonalize  
%     [U,~] = qr(AJ,0);
%     [V,~] = qr(AI',0);
    [U,~,~] = svd(AJ, 'econ');
    [V,~,~] = svd(AI', 'econ');
    
    
    S_new = 0;
    old_rank = 0;
    
    for iter = 2:Max_iter
        
        [J_hat, rank_hat] = sample_cols(V, J, rank, [], maxvol_tol);
        I_hat             = sample_rows(I, dims, rank_hat, n0, iter, idxS);
%         I_hat             = sample_rows_ql(I, dims, rank_hat, n0, iter, idxS);
%        I_hat = sample_rows_maxvol(U, I, rank);
        
        if rank_hat == 0
            break;
        end
        
        I = union(I, I_hat);
        J = union(J, J_hat);
        
        % sample selected columns and rows       
        U_aux = fU(J_hat);
        V_aux = fV(I_hat);
        
%         U_aux = U_aux - U * (U' * U_aux);
%         
%         [U_aux,~,~] = svd(U_aux,'econ'); 
        
        % form a new submatrix
        U = [U, U_aux];
        V = [V, V_aux'];
                
%         [U,~] = qr(U,0);
%         [V,~] = qr(V,0);
        [U,~,~] = svd(U, 'econ');
        [V,~,~] = svd(V, 'econ');


        U_hat = U(I,:);
        V_hat = V(J,:);
        A_hat = fA_hat(I, J);
        
        [Qu_hat,Su_hat,Vu_hat] = truncated_svd(U_hat, 1e-15); 
        [Qv_hat,Sv_hat,Vv_hat] = truncated_svd(V_hat', 1e-15);
        
        M = Su_hat \ (Qu_hat' * A_hat * Vv_hat) / Sv_hat;
        
        [Qm,Sm,Vm] = truncated_svd(M, tol);        
        
        Qu = Vu_hat * Qm;
        Qv = Sm * Vm' * Qv_hat';
        
        S_old = S_new; 
        S_new = diag(Sm);
        rank = find_rank_from_S(S_new, tol);
        
        max_rank = max(max_rank,rank);
        
        if abs(rank - old_rank) < 2 || ((rank - old_rank)< 0)
            stall = stall + 1;
        end
        

        
%         U = U(:,1:max_rank);
%         V = V(:,1:max_rank);
              
        % compute norm of difference in singular values (old - new)
        if f_term_criteria(S_old, S_new, old_rank, rank, tol) || (stall == 10)
            X = U * Qu;
            W = Qv * V';
            return;
        end
                
       old_rank = rank;
       
      
       
    end
end