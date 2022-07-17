function I_hat = sample_rows(I, dims, rank, n0, iter, idxS)


n = n0 + 2*(iter - 1);

x_cheb = @(n) real(exp(1i * linspace(0,pi,n + 1)));

ix = unique(nonzeros(ceil(dims.L_vie/2 * (1 + x_cheb(n)))));
iy = unique(nonzeros(ceil(dims.M_vie/2 * (1 + x_cheb(n)))));
iz = unique(nonzeros(ceil(dims.N_vie/2 * (1 + x_cheb(n)))));

[jY,iX,kZ] = meshgrid(iy, ix, iz);

% dom_dims = dims.L_vie * dims.M_vie * dims.N_vie;
dom_dims = dims.N_scat;

idx = sub2ind(dims.vie(1:3), iX(:), jY(:), kZ(:));

idx = intersect(idx, idxS);
idx = find(ismember(idxS, idx));

I_hat = repmat(idx, dims.ql,1) + kron(dom_dims * (0:dims.ql-1).', ones(size(idx)));
I_hat = setdiff(I_hat, I);

sz_I_hat = size(I_hat,1);
i_loc = sort(randsample(sz_I_hat, min(rank, sz_I_hat)));
I_hat = I_hat(i_loc);
