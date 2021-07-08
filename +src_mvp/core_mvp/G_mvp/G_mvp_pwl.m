function Jout = G_mvp_pwl(Jin, dx, dom_dims, basis_dims)

% reshape input currents to tensor; allocate memory for Jout
Jin_tensor  = reshape(Jin,dom_dims);

% allocate memory for output results
Jout = zeros(dom_dims, 'like', Jin_tensor);

% voxel's volume
V = dx^3;

% Gram matrix for single voxel
G_vox = [ V   0     0     0
          0   V/12  0     0
          0   0     V/12  0
          0   0     0     V/12];

% vectorial Gram matrix
G = kron(eye(basis_dims.q), G_vox);

%% matrix vector product with G
for i = 1:dom_dims(end)
    Jout(:,:,:,i) = G(i,i) * Jin_tensor(:,:,:,i);
end

