function Jout =  N_mvp_pwl_cell(Jin, fN_cell, N_cell_idx, dom_dims, op_dims, dims)

% get dimensions of a domain
L = dom_dims(1);
M = dom_dims(2);
N = dom_dims(3);

% reshape/allocate currents
Jin_tensor = reshape(Jin, dom_dims);
Jout       = Jin_tensor;
Jin_fft    = zeros([op_dims(1:3), dom_dims(end)], 'like', Jin);
N_Jin_ifft = Jin_fft;
    
% compute fftn's for current tensors
for i=1:dom_dims(end)
    Jin_fft(:,:,:,i) = fftn(Jin_tensor(:,:,:,i), op_dims(1:3));
end

% sgn_vec  = N_cell_idx.sgn_vec;
left_vec = N_cell_idx.left_vec; 

% tic;
%% compute fields produced by PWL basis
for i=1:dims.ql

    % get appropriate operator components
%     i_sgn = sgn_vec(:,i);
    
    % compute element-wise product between fft of N operator and a current
%     N_Jin = reshape(reshape(fN_cell{i} .* Jin_fft,[],12) * i_sgn, op_dims(1:3));
    N_Jin = sum(fN_cell{i} .* Jin_fft,4);
    
    % get resulting fields
    N_Jin_ifft(:,:,:,i) = ifftn(N_Jin);
end

% toc;
%% original implementation

% get results
Jout(:,:,:,left_vec) = N_Jin_ifft(1:L,1:M,1:N,1:dims.ql);

% reshape to Jin dimensions 
Jout = reshape(Jout, size(Jin));

