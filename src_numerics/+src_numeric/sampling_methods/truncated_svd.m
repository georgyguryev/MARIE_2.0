function [U,S,V] = truncated_svd(A, trunc_tol)

[U,S,V] = svd(A, 'econ');

s_vec = diag(S);
        
rank = find_rank_from_S(s_vec, trunc_tol);

U = U(:,1:rank);
S = S(1:rank, 1:rank);
V = V(:,1:rank);
       