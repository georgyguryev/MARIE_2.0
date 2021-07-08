function [PS] = projection_assembly_vox(fP_mvp, coil, scatterer, projector, dims, task_settings, freq)
% function [P] = projection_assembly(,,,)  - assembles the projection
% matrix that maps surface currents to volumetric fictitious currents;
% this implementation relies on voxel integration as opposed to simple
% point-matching, implemented in project_RWG_on_PWX();
% ----------------------------------------------------------------------- %
% INPUT parameters: 
%       - fN_mvp_near(Jin) MVP handle to apply N operator to 'near' domain
%       - COIL struct. with RWG basis infornation
%       - RHBM struct. with geometric and material properties
%       - PROJ struct. with extended domain info
%       - freq - working frequency
% ----------------------------------------------------------------------- %
% OUTPUT:
%       P - pojection matrix
% ----------------------------------------------------------------------- %

%% setup constants/parameters

% electromagnetic constants and frequency
emu = src_utils.EM_utils(freq);

% copy centers of grid cubes
r             = scatterer.dom_ext.r;
res           = scatterer.dom_ext.res;
idxS_ext      = scatterer.index_ext.S_1d;

% get projection-related parameters
frame_width = projector.near_boundary_width;
Cell_Ctr    = projector.central_cell_id;
Cells       = projector.cell_id;

N_sie      = dims.N_sie;
Next       = dims.Nvox_ext;
N_exp_3D   = dims.N_exp_3D;
N_near_1D  = dims.N_near_1D;
N_near_3D  = dims.N_near_3D;

% select appropriate basis size
PWX_sz  = dims.ql;

N_ext_dofs = Next * dims.ql;
N_exp_dofs = N_exp_3D * dims.ql;
N_vie_dofs = dims.N_scat * dims.ql;
N_collocation_dofs = (N_near_3D - (N_near_1D - 2 * frame_width).^3) * dims.ql;

% resize resulting projection matrix with expansion coefficients
P = sparse(N_ext_dofs, N_sie);

% for collocation at the boundary voxels
in_near_st  = frame_width + 1;
in_near_end = N_near_1D - frame_width;

% get coordinates of a 'near zone' domain
r_near    = r(1:N_near_1D, 1:N_near_1D, 1:N_near_1D,:);
r_near_in = r(in_near_st:in_near_end,...
              in_near_st:in_near_end,...
              in_near_st:in_near_end, :);


% reshape r_near and r_near_in to extract boundary coordinates
r_near    = reshape(r_near, N_near_1D^3, dims.q);
r_near_in = reshape(r_near_in, (N_near_1D - 2)^3, dims.q);
% r_near_in = reshape(r_near_in, (N_near_1d - 2 * frame_width)^3, PWX_sz);

% get expansion dofs indices
[dofs_idx, I_Jb_idx, J_Jb_idx] = src_utils.get_exp_dofs(Cells,idxS_ext.',PWX_sz, Next);

V_Jb_idx = ones(size(I_Jb_idx));

S = sparse(I_Jb_idx, J_Jb_idx, V_Jb_idx, N_ext_dofs, N_vie_dofs);

%% Check order !!!
r_near_bound = setdiff(r_near, r_near_in, 'row', 'stable');

max_I = size(r_near_bound,1);
ind_r_near = randi(max_I, N_exp_3D, 1);
r_near_bound = r_near_bound(ind_r_near,:);

N_collocation_dofs = max_I;

% r_near_bound = squeeze(r(N_near_1d,N_near_1d,N_near_1d,:)).';

% set a displacement vector (displacement of central vox with respect to
% the left bottom voxel
N_disp = floor(N_near_1D / 2); 
r_disp = - squeeze(r(1,1,1,:)).' - res * [N_disp, N_disp, N_disp];

%% assemble matrix A

A = zeros(N_exp_dofs, N_exp_dofs);

% form an explicit matrix representation
for i_exp = 1: N_exp_dofs
    e_i = zeros(N_exp_dofs,1);
    e_i(i_exp) = 1;
    E_vox_bound = fP_mvp(e_i);
    A(1:N_exp_3D,i_exp) = E_vox_bound(ind_r_near);
    A(N_exp_3D + 1:2 * N_exp_3D, i_exp) = E_vox_bound(ind_r_near + max_I);
    A(2 * N_exp_3D + 1:3 * N_exp_3D, i_exp) = E_vox_bound(ind_r_near + 2 * max_I);

end


% multiply by 1/ce scaling factor
A = 1 / emu.ce * A;

%% set up quadratures 

% get order of quadratures
Quad_order_sie  = task_settings.vsie.Np_quad_coup_sie;
Quad_order_vie  = task_settings.vsie.Np_quad_coup_vie;

% get quadrature points and weights
[wp_vie, z_vie]                          = gauss_1d(Quad_order_vie);
[Np_sie, z1_sie, z2_sie, z3_sie, wp_sie] = dunavant_rule(Quad_order_sie);

% form a vector of quadratures
VIE_quads = [Quad_order_vie; wp_vie; z_vie];
SIE_quads = [Np_sie; wp_sie; z1_sie; z2_sie; z3_sie];

%% select radius vectors for boundary voxels for each

Z_bc = zeros(N_exp_dofs, N_sie);

% profile on;

for i_sie = 1:N_sie
    
     % get index of current central expansion voxel
    ctr_vox_idx = Cell_Ctr(:, i_sie);
%     cel_vox_idx = Cells(i_sie, :);    % returns array of global indices of cell voxels
    
    % get radius vector of central voxel
    r_c = (squeeze(r(ctr_vox_idx(1), ctr_vox_idx(2), ctr_vox_idx(3),:))).';
    
    % get total displacement of boundary domain 
    r_full_disp = r_disp + r_c;
    
    r_bound     = r_near_bound + r_full_disp;
       
    % coumpute right-hand side/ excitation due to current RWG basis
%     Z_bc(:,i_sie) = Assembly_COUP_rwg(r_bound,coil,dims,task_settings,...
%                                       freq,res,i_sie); 

    r_rwg = src_coupling.get_rwg_vertices(coil, i_sie);
    
        
    Z_bc(:,i_sie) = Assemble_rwg_coupling(r_bound.', r_rwg, SIE_quads, VIE_quads, res, dims.l, emu.k0);
                                  
end
% profile off;
% profile viewer;


rhs = res.^3 * Z_bc;
    
Jsol = A \ rhs;


%% fill projection matrix
tic
for i_sie = 1:N_sie
    
    % get current dofs indicies
    dofs_sie_idx = dofs_idx(:, i_sie);
    
    % 
    P(dofs_sie_idx,i_sie) = Jsol(:,i_sie);
    
end
toc;

keyboard;
% concatinate projector matrices
PS = [P,S];

% keyboard;









