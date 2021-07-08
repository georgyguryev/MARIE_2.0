function [idx,prop_ext] = get_scat2basis_indices(basis_dom, basis_ind, scatterer, freq)
% function get_scat2basis_indices() maps scatterer voxel indices to the
% ones that correspond to the same voxels within the basis body model

% get scatterer and basis body domains
dom_vie     = scatterer.dom_vie;
dom_basis   = basis_dom;

index_basis = basis_ind;

% allocate indices for material properties
prop_ext.epsilon_r = ones(size(dom_basis.x_tensor));
prop_ext.sigma_e   = zeros(size(dom_basis.x_tensor));
prop_ext.rho       = zeros(size(dom_basis.x_tensor));

% mapping indexes from extended to VIE domain
idx_vie_x = find((dom_basis.x - min(dom_vie.x))>-eps& ...
                 (dom_basis.x - max(dom_vie.x))<eps);
idx_vie_y = find((dom_basis.y - min(dom_vie.y))>-eps& ...
                 (dom_basis.y - max(dom_vie.y))<eps);
idx_vie_z = find((dom_basis.z - min(dom_vie.z))>-eps& ...
                 (dom_basis.z - max(dom_vie.z))<eps);
                          
% map material properties to the extended grid
prop_ext.rho(idx_vie_x,idx_vie_y,idx_vie_z)        = scatterer.prop_vie.rho(:,:,:);
prop_ext.sigma_e (idx_vie_x,idx_vie_y,idx_vie_z)   = scatterer.prop_vie.sigma_e(:,:,:);
prop_ext.epsilon_r (idx_vie_x,idx_vie_y,idx_vie_z) = scatterer.prop_vie.epsilon_r(:,:,:);
 
%%

% compute omega
emu = src_utils.EM_utils(freq);

% define complex relative permittivity and susceptibility
Mr  = prop_ext.epsilon_r + prop_ext.sigma_e / (emu.ce);

% get 1-D indices that correspond to voxels 
idx_S_1D = find(abs(Mr - 1) > 1e-12);

idx_basis = index_basis.S_1d;

idx = ismember(idx_basis, idx_S_1D);

if nnz(idx) < size(idx_S_1D,1)
    error("Current body is fully or partially outside of the basis bounding box! ");
end

