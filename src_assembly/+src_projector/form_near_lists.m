function [C2Exp_centers, RWG_to_PWX_cells, C2Vox_near_lists,...
          C2C_near_list, C2B_near_list, index_exp2near] = form_near_lists(coil, scatterer, dims)    
% function looks for the closest voxel expansion cell 
% for a given RWG basis
%  [RWG_to_PWX_centers, RWG_to_PWX_cells, RWG_PWX_list_near, RWG_RWG_list_near]

% get dimensions of coil and expansion domains
N_sie    = dims.N_sie;
N_exp_3D = dims.N_exp_3D;


% extract coordinates of (terminal + inner) RWG centers
[RWG_cntr] = get_RWG_centers(coil);

%% 

% allocate memory
C2Exp_centers    = zeros(3, N_sie);
C2C_near_list    = containers.Map;
C2B_near_list    = containers.Map;

C2Vox_near_lists = struct('ext', zeros(dims.N_near_3D, N_sie),...
                          'near', {cell(N_sie,1)});


%% find nearest center of expansion voxel (for given edge)
for i_sie = 1:N_sie
    
    % get current RWG center
    r_sie = RWG_cntr(i_sie,:).';
    
    % associate current RWG with the nearest voxel
    C2Exp_centers(:,i_sie)  = src_projector.find_nearest_voxel(r_sie, scatterer);
    
end
 
% assign an expansion cell for given RWG function
% RWG_to_PWX_cells(:,i_sie) = src_projector.assign_expansion_cell(C2Exp_centers, dims);
RWG_to_PWX_cells = src_projector.assign_expansion_cell(C2Exp_centers, dims, 'ext');

% form a simple near interaction list with indexing within extended domain
[C2Vox_near_lists.ext, C2Vox_near_lists.near] = src_projector.form_expansion_near_list(C2Exp_centers, dims);

% check for overlap between scatterer voxels and near domain
I_sn = ismember(C2Vox_near_lists.ext(:), scatterer.index_ext.S_1d);

% reshape mask to matrix
I_sn = reshape(I_sn, dims.N_near_3D, N_sie);

for i_sie = 1:dims.N_sie
% if current RWG basis has scatterer voxel in near zone - add to list
    I_mask = I_sn(:,i_sie);
    if nnz(I_mask) > 0
        C2B_near_list(num2str(i_sie)) = C2Vox_near_lists.ext(I_mask,i_sie);
    end
end


%% form C2C near list

a = C2Exp_centers;

for i_sie = 1:N_sie
    
    dist_3D  = abs(bsxfun(@minus, a, a(:,i_sie)));
    dist_inf =  max(dist_3D);
    
    % find near C2C interacting pairs
    near_list = find(dist_inf < ((dims.N_near_1D - 1)/2));
    
    % update the near list 
    C2C_near_list(num2str(i_sie)) = near_list;
        
end

%% define mapping from global index of expansion voxels to local near index

N_start = (dims.N_near_1D - dims.N_exp_1D) / 2;

% get indices of voxels
exp_vox_1D = N_start + (1:dims.N_exp_1D);

% generate indices for all expansion voxels within near zone
[i_exp_L,j_exp_M,n_exp_N] = meshgrid(exp_vox_1D, exp_vox_1D, exp_vox_1D);


index_1D = reshape(sub2ind(dims.near(1:3), j_exp_M, i_exp_L, n_exp_N), N_exp_3D,1);

% form local index for all ql components of expansion basis in near basis
index_exp2near = src_scatterer.Scatterer_index(index_1D, dims.N_near_3D);


end