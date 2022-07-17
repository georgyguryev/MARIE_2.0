function [I_exp_vox] = assign_expansion_cell(r_exp_ctr, dims, indexing_domain)

I_exp_vox = zeros(dims.N_exp_1D.^3, dims.N_sie);

N_exp_1D = dims.N_exp_1D;

switch indexing_domain
    case 'ext'
        L = dims.L_ext;
        M = dims.M_ext;
        N = dims.N_ext;
    case 'near'
        L = dims.N_near_1D;
        M = L;
        N = L;
end

% we assume that Ncell is n_1d^3 , where n_1d is number of cubes in 1 dim
idx_range = -(N_exp_1D - 1) / 2:1:(N_exp_1D-1) / 2;
idx_range = repmat(idx_range, dims.N_sie,1);

idx_x = repmat(r_exp_ctr(1,:).', 1, N_exp_1D);
idx_y = repmat(r_exp_ctr(2,:).', 1, N_exp_1D);
idx_z = repmat(r_exp_ctr(3,:).', 1, N_exp_1D);

% form the range of voxel indices
idx_cell_x = (idx_range + idx_x).';
idx_cell_y = (idx_range + idx_y).';
idx_cell_z = (idx_range + idx_z).';

for i_sie = 1:dims.N_sie
    
    [X,Y,Z] = meshgrid(idx_cell_x(:,i_sie), idx_cell_y(:,i_sie), idx_cell_z(:,i_sie));
    
    % map 3D indices to 1D
    exp_glob_index = sub2ind([L M N], X, Y, Z);
    
    exp_glob_index = permute(exp_glob_index, [2 1 3]);
    
    I_exp_vox(:,i_sie) = exp_glob_index(:);

end
