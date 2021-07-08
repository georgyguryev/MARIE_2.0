function [I_near_glob, I_near_loc] = form_expansion_near_list(r_exp_ctr, dims)
% function [idx_near_list] = form_expansion_near_list(idx_x, idx_y, idx_z, Ncell, r)
% function forms a list of "near" interactions between PWX expansion voxels


% we assume that Ncell is n_1d_near^3 , where n_1d_near is number of cubes
% per Cell; for now the total number of near expansion voxels 1D is given by 
N_near_1D  = dims.N_near_1D;

% domain dimensions 
L = dims.L_ext;
M = dims.M_ext;
N = dims.N_ext;

I_near_glob = zeros(dims.N_near_3D, dims.N_sie);
I_near_loc  = cell(dims.N_sie,1);

% we assume that Ncell is n_1d^3 , where n_1d is number of cubes in 1 dim
idx_range = -ceil((N_near_1D-1)/2):ceil((N_near_1D-1)/2);
idx_range = repmat(idx_range, dims.N_sie,1);

idx_x = repmat(r_exp_ctr(1,:).', 1, N_near_1D);
idx_y = repmat(r_exp_ctr(2,:).', 1, N_near_1D);
idx_z = repmat(r_exp_ctr(3,:).', 1, N_near_1D);


% form the range of voxel indices
idx_cell_x = (idx_range + idx_x).';
idx_cell_y = (idx_range + idx_y).';
idx_cell_z = (idx_range + idx_z).';


% mask indices that are out of range
idx_cell_x(idx_cell_x >  L) = NaN;
idx_cell_y(idx_cell_y >  M) = NaN;
idx_cell_z(idx_cell_z >  N) = NaN;

idx_cell_x(idx_cell_x < 1) = NaN;
idx_cell_y(idx_cell_y < 1) = NaN;
idx_cell_z(idx_cell_z < 1) = NaN;

% % generate mesh for current near domain 
% [X,Y,Z] = meshgrid(idx_cell_x, idx_cell_y, idx_cell_z);
% 
% % map 3D indices to 1D
% glob = sub2ind([L M N], X, Y, Z);
% 
% % permute indices for appropriate dimension
% glob = permute(glob, [2 1 3]);
% 
% % replace NaN's with zeros (indicates that some voxels are outside ext domain)
% glob(isnan(glob)) = 0;
% 
% % return global indices
% idx_near_list = glob(:).';


for i_sie = 1:dims.N_sie
    
    [X,Y,Z] = meshgrid(idx_cell_x(:,i_sie), idx_cell_y(:,i_sie), idx_cell_z(:,i_sie));
    
    % map 3D indices to 1D
    exp_glob_index = sub2ind([L M N], X, Y, Z);
    
    exp_glob_index = permute(exp_glob_index, [2 1 3]);
    
    % replace NaN's with zeros (indicates that some voxels are outside ext domain)
    exp_glob_index(isnan(exp_glob_index)) = 0;
    
    I_near_glob(:,i_sie) = exp_glob_index(:);
    I_near_loc{i_sie} = find(exp_glob_index(:));

end

