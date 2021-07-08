function P_j = fill_P_par(Ic_new, k_vie, j_vie, N_cols)

% fill with zeros
P_j = sparse(1,N_cols);

% fill projection entries 
P_j(1,j_vie) = Ic_new(k_vie);

