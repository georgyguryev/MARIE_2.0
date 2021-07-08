function [U_N, U_K] = generate_basis_CS(task_settings, dims, scatterer, coil, freq)

N_samples = 150; 
tol = task_settings.basis.Tolerance;

Ct = coil.Ct;
P  = zeros(size(Ct,2),1);

idxS = scatterer.index_vie.S_3d;

scat_r = reshape(scatterer.dom_vie.r(idxS), dims.N_scat,dims.q);


delta = @(r) bsxfun(@minus, scat_r.', r);
dist = @(d) d(1,:).^2 + d(2,:).^2 + d(3,:).^2;

for i = 1:dims.N_sie
    
    i_index = find(coil.index == i);
    i_nodes  = coil.edge(:,i_index);
    r1 = coil.node(:,i_nodes(1));
    r2 = coil.node(:,i_nodes(2));
    
    delta_1 = delta(r1);
    delta_2 = delta(r2);
    ave_dist = dist(delta_1) + dist(delta_2);
    P(i) = 1 / min(ave_dist(:));
end


cols = 1:dims.N_sie;

P = P / sum(P);
m_1 = find(mnrnd(N_samples, P));
% m_1 = randi(length(P), N_samples, 1);

% samp_cols = datasample(cols,N_samples,'Replace',false);
samp_cols = cols(m_1);
% cols = setdiff(cols,samp_cols);

[Zbc_N_sub, Zbc_K_sub] = src_coupling.assemble_coupling_matrices(task_settings, dims,...
                                                   scatterer, coil, freq, samp_cols);                                               
                                               
P(m_1)    = 0;
P = nonzeros(P) / sum(nonzeros(P));                                               
                                               
[U_N_sub,S,~] = svd(Zbc_N_sub,'econ');
[U_K_sub,~,~] = svd(Zbc_K_sub,'econ');

Z_N = Zbc_N_sub;
Z_K = Zbc_K_sub; 
clear Zbc_N_sub Zbc_K_sub;


s_new = diag(S);
s_half = sum(s_new(ceil(end/2)+1:end) / s_new(1));

U_N = U_N_sub;
U_K = U_K_sub;
clear U_N_sub U_K_sub;

l = ceil(size(U_N,2));
gamma = sqrt(dims.N_sie / l);

while (size(U_N,2) < dims.N_sie) && (gamma * s_half > tol)
    
    m_1 = find(mnrnd(N_samples, P));
    
%     m_1 = randi(length(P), N_samples, 1);

%     samp_cols = datasample(cols,N_samples,'Replace',false);
    samp_cols = cols(m_1);
%     cols = setdiff(cols,samp_cols);
    
    [Zbc_N_sub, Zbc_K_sub] = src_coupling.assemble_coupling_matrices(task_settings, dims,...
                                                                   scatterer, coil, freq,...
                                                                   samp_cols);

    Z_N = [Z_N Zbc_N_sub];
    Z_K = [Z_K Zbc_K_sub];
    
    clear Zbc_N_sub Zbc_K_sub;
    
    [U_N,S,~] = svd(Z_N, 'econ');
    
    P(m_1)    = 0;
    P = nonzeros(P) / sum(nonzeros(P));
    
    
    clear U_sub;
    
    s_new  = diag(S);
    s_half = sum(s_new(ceil(end/2)+1:end) / s_new(1));
    
    l = ceil(size(U_N,2));
    
    gamma = sqrt(dims.N_sie / l);
        
end

[U_K,~,~] = svd(Z_K, 'econ');


U_N = U_N(:,1:ceil(end/2));
U_K = U_K(:,1:ceil(end/2));

