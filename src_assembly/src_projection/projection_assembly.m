function [P,S] = projection_assembly(fP_mvp, COIL, RHBM, PROJ, freq)
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
EM_params;

% get coil-related data
index = COIL.index;
etod  = COIL.etod;
node  = COIL.node;
elem  = COIL.elem;


% get gemetric information from RHBM
r             = RHBM.r;
idxS_ext      = RHBM.idxS_ext.';

% get projection-related parameters
PWX_flag      = PROJ.PWX_flag;
GL_order_sie  = PROJ.GL_order_sie;
GL_order_vie  = PROJ.GL_order_vie;
constr_order  = PROJ.constr_order;
N_col_pts     = PROJ.N_col_pts;
col_cube_side = PROJ.Col_dist;      % vector of collocation distances
Cell_Ctr      = PROJ.Cell_Centers;
Cells         = PROJ.Cells;
N_vox_1D      = PROJ.Nexp_1D;
Near_PWX      = PROJ.Near_PWX; 
frame_width   = PROJ.frame_width;


% get dimensions of extended VSIE domain
[Lx, My, Nz, ~] = size(r);


N_sie     = max(COIL.index);
N_ext     = Lx * My * Nz;
N_exp     = size(Cells,2);
N_near    = size(Near_PWX,2);
N_near_1d = round((N_near)^(1/3));

% select appropriate basis size
PWX_sz       = pwx_size(PWX_flag);
N_exp_dofs   = N_exp * PWX_sz;
N_ext_dofs   = N_ext * PWX_sz;
N_near_dofs  = N_near * PWX_sz;
N_near_bound = N_near_dofs - PWX_sz * (N_near_1d - 2 * frame_width)^3;
N_colloc_sz  = N_near_dofs - PWX_sz * (N_near_1d - 2)^3;
N_vie_dofs   = size(idxS_ext,2) * PWX_sz;
% N_colloc_sz = PWX_sz;

% resize resulting projection matrix with expansion coefficients
P = sparse(N_ext_dofs, N_sie);
S = sparse(N_ext_dofs, N_vie_dofs);

dofs_idx = zeros(N_exp_dofs,N_sie);

x = r(:,:,:,1);
y = r(:,:,:,2);
z = r(:,:,:,3);

res = x(2,1,1) - x(1,1,1);

% % for collocation in the frame voxels
% in_near_st  = frame_width + 1;
% in_near_end = N_near_1d - frame_width;

% for collocation at the boundary voxels
in_near_st  = 2;
in_near_end = N_near_1d - 1;

% get coordinates of a 'near zone' domain
r_near    = r(1:N_near_1d, 1:N_near_1d, 1:N_near_1d,:);
r_near_in = r(in_near_st:in_near_end,...
              in_near_st:in_near_end,...
              in_near_st:in_near_end, :);


% reshape r_near and r_near_in to extract boundary coordinates
r_near    = reshape(r_near, N_near_1d^3, pwx_size(0));
r_near_in = reshape(r_near_in, (N_near_1d - 2)^3, pwx_size(0));
% r_near_in = reshape(r_near_in, (N_near_1d - 2 * frame_width)^3, PWX_sz);

% get expansion dofs indices
[dofs_idx, I_Jb_idx, J_Jb_idx] = get_exp_dofs(Cells,idxS_ext,PWX_sz, N_ext);

V_Jb_idx = ones(size(I_Jb_idx));

S = sparse(I_Jb_idx, J_Jb_idx, V_Jb_idx, N_ext_dofs, N_vie_dofs);

%% Check order !!!
r_near_bound = setdiff(r_near, r_near_in, 'row', 'stable');

% r_near_bound = squeeze(r(N_near_1d,N_near_1d,N_near_1d,:)).';

% set a displacement vector (displacement of central vox with respect to
% the left bottom voxel
N_disp = floor(N_near_1d / 2); 
r_disp = - squeeze(r(1,1,1,:)).' - res * [N_disp, N_disp, N_disp];

%% assemble matrix A

A = zeros(N_colloc_sz, N_exp_dofs);

% form an explicit matrix representation
for i_exp = 1: N_exp_dofs
    e_i = zeros(N_exp_dofs,1);
    e_i(i_exp) = 1;
    A(:,i_exp) = fP_mvp(e_i);
end


% multiply by 1/ce scaling factor
A = 1 / ce * A;

% Jin = zeros(N_near_dofs,1);
% Jin(4 + 4 * 11 + 4 * 121) = 1;
% 
% rhs = 1 / ce * fP_mvp_test(Jin);
% 
% Jsol = A \ rhs;

tic 
%% select radius vectors for boundary voxels for each 

Z_bc = zeros(N_colloc_sz, N_sie);

for i_sie = 1:N_sie
    
     % get index of current central expansion voxel
    ctr_vox_idx = Cell_Ctr(i_sie, :);
%     cel_vox_idx = Cells(i_sie, :);    % returns array of global indices of cell voxels
    
    % get radius vector of central voxel
    r_c = (squeeze(r(ctr_vox_idx(1), ctr_vox_idx(2), ctr_vox_idx(3),:))).';
    
    % get total displacement of boundary domain 
    r_full_disp = r_disp + r_c;
    
    r_bound = r_near_bound + r_full_disp;
       
    % coumpute right-hand side/ excitation due to current RWG basis
    
    %%  TODO!!! Adjust assembly function to compute a vector Zbc for given i_sie; remove 1/ce normalization
%     z_bc = Assembly_SCOUP_QMEX_par(r_bound,index,etod,node,elem,freq,GL_order_sie);
    

    Z_bc(:,i_sie) = Assembly_COUP_rwg(r_bound,index,etod,node,elem,PWX_flag,...
                                      GL_order_sie,GL_order_vie,freq,res,i_sie);    
end

rhs = res^3 * Z_bc;
    
Jsol = A \ rhs;


%% fill projection matrix

for i_sie = 1:N_sie
    
    % get current dofs indicies
    dofs_sie_idx = dofs_idx(i_sie,:);
    
    % 
    P(dofs_sie_idx,i_sie) = Jsol(:,i_sie);
    
end

toc

% keyboard;









