function Jout =  N_mvp_pwl(Jin, fN, dom_dims, op_dims, dims)

% get dimensions of a domain
L = dom_dims(1);
M = dom_dims(2);
N = dom_dims(3);

% reshape/allocate currents
Jin_tensor = reshape(Jin, dom_dims);
Jout       = Jin_tensor;
Jin_fft    = zeros([op_dims(1:3), dom_dims(end)], 'like', Jin);


% mappings that exploit symmetry in diadyc Green's Funciton and interacting
% faces of voxels
map_pq   = [1 2 3; 2 4 5; 3 5 6];
map_lplq = [1 5 6 7; -5 2 8 9; -6 8 3 10; -7 9 10 4];
  
% compute fftn's for current tensors
for i=1:dom_dims(end)
    Jin_fft(:,:,:,i) = fftn(Jin_tensor(:,:,:,i), op_dims(1:3));
end

% generate indices to unfold loops
[p_msh, lp_msh] = meshgrid(1:dims.q, 1:dims.l);

p_vec  = kron(p_msh(:), ones(dims.ql,1));
lp_vec = kron(lp_msh(:), ones(dims.ql,1));
q_vec  = repmat(p_msh(:), dims.ql,1);
lq_vec = repmat(lp_msh(:), dims.ql,1);
pq_vec     = sub2ind([dims.q, dims.q], p_vec, q_vec);
lplq_vec   = sub2ind([dims.l, dims.l], lp_vec, lq_vec);
left_vec   = sub2ind([dims.l, dims.q], lp_msh, p_msh);

% define final indexing for N_mvp_pwl
idx_pq_vec   = reshape(map_pq(pq_vec), [dims.ql, dims.ql]);
idx_lplq_vec = reshape(map_lplq(lplq_vec), [dims.ql, dims.ql]);
idx_lq_vec   = reshape(sub2ind([10, 6], abs(idx_lplq_vec), idx_pq_vec), dims.ql, dims.ql);
sgn_vec      = sign(idx_lplq_vec);

% allocate memory for temporary results (to be cropped from FFT domain to domain dims)
N_Jin_ifft = Jin_fft;

%% compute fields produced by PWL basis
for i=1:dims.ql

    % get appropriate operator components
    i_lq = idx_lq_vec(:,i);
    
    Nij = fN(:,:,:,i_lq);
    i_sgn = sgn_vec(:,i);
        
    % compute element-wise product between fft of N operator and a current
    N_Jin = reshape(reshape(fN(:,:,:,i_lq) .* Jin_fft,[],12) * i_sgn, op_dims(1:3));
     
%     Nij = fN(:,:,:,i_lq);
                
%     N_Jin = reshape(reshape(Nij .* Jin_fft,[],12) * i_sgn, op_dims(1:3));
    
    % get resulting fields
    N_Jin_ifft(:,:,:,i) = ifftn(N_Jin);

end

%% original implementation

% get results
Jout(:,:,:,left_vec) = N_Jin_ifft(1:L,1:M,1:N,1:dims.ql);

% reshape to Jin dimensions 
Jout = reshape(Jout, size(Jin));


