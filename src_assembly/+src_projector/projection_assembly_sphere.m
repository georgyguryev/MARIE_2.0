function [PS] = projection_assembly_sphere(coil, scatterer, projector, dims, task_settings, freq)
% function [] = project_RWG_on_PWX(COIL, RHBM, PROJ, f)
% function finds expansion coefficients ("projection") of RWG basis
% on PWX expansion cell

% 16.01.18 Commented part for projection with linear constraints
%          added collocation points at different distance



% 25.04.2021 Updated code to assemble Electric and Magnetic projection
% matrices

%% setup constants/parameters

% get number of concentric collocation spheres

% copy centers of grid cubes
r             = scatterer.dom_ext.r;
res           = scatterer.dom_ext.res;
idxS_ext      = scatterer.index_ext.S_1d;
col_distance  = (dims.N_near_1D -1) / 2 + 1;
Cell_Ctr      = projector.central_cell_id;
Cells         = projector.cell_id;

Nexp_1D    = dims.N_exp_1D;
N_sie      = dims.N_sie;
Next       = dims.Nvox_ext;
N_cells    = dims.N_exp_3D;



%% select number of quadrature points based on dimension of the expansion domain

switch dims.N_exp_1D
    % select number of quad points on a sphere
    case 3
        N_colloc_pts = 86;
    case 5
        N_colloc_pts = 86;
    case 7
        N_colloc_pts = 230;
    otherwise
        N_colloc_pts = 230;
end

% get coil data
index = coil.index;
etod  = coil.etod;
node  = coil.node;
elem  = coil.elem;

% get electromagnetic constants
emu = src_utils.EM_utils(freq);

% select appropriate basis size
PWX_sz  = dims.ql;
N_basis = dims.l;

Next_dofs  = Next * dims.ql;
N_vie_dofs = dims.N_scat * dims.ql;


%% set up quadratures 

% get order of quadratures

Quad_order_C2Col = 1;
Quad_order_sie   = task_settings.vsie.Np_quad_coup_sie;
Quad_order_vie   = task_settings.vsie.Np_quad_coup_vie;

% get quadrature points and weights
[wp_vie, z_vie]                          = gauss_1d(Quad_order_C2Col);
[Np_sie, z1_sie, z2_sie, z3_sie, wp_sie] = dunavant_rule(Quad_order_sie);

% form a vector of quadratures
VIE_quads = [Quad_order_C2Col; wp_vie; z_vie];
SIE_quads = [Np_sie; wp_sie; z1_sie; z2_sie; z3_sie];

%% allocate memory for projection values and indices
Val   = cell(N_sie,1);
I_val = cell(N_sie,1);

%% Assemble cell projector martix


% tic;

r_o = [0; 0; 0]; 

% generate collocation points 
col_points = src_assembly.generate_collocation_points(res, col_distance, N_colloc_pts);

% generete cube centers for elementary Cell
[vox_centers] = src_assembly.generate_cube_centers(r_o, res, Nexp_1D);

% compute Lop/Nop matrix that maps volumetrix currents to E in colloc. points      
% [A,B] = E_DGF_PWX_Col(vox_centers.', res, col_points, Quad_order_vie, freq, dims);
[A,B] = E_DGF_PWX_Col(vox_centers.', res, col_points, 3, freq, dims);


% allocate memory for the rhs
rhs_E = zeros(3 * N_colloc_pts, N_sie);
rhs_H = zeros(3 * N_colloc_pts, N_sie);

for i_sie = 1:N_sie

    % get index of current central expansion voxel
    ctr_vox_idx = Cell_Ctr(:, i_sie);
    
    % get coordinates of central expansion voxels and get col. pts
    r_c = squeeze(r(ctr_vox_idx(1), ctr_vox_idx(2), ctr_vox_idx(3),:));
    
    % form current collocation points
    cur_col_points = col_points + r_c;  % + repmat(r_c, 1, size(col_points,2));
        
    % get coordinates of current rwg    
    rp_i = coil.rp(:,i_sie);
    rn_i = coil.rn(:,i_sie);
    r2_i = coil.r2(:,i_sie);
    r3_i = coil.r3(:,i_sie);
    
    % get fields produced by unitary excitation of current rwg triangle
    [rhs_E_i, rhs_H_i] = Assemble_rwg_coupling_matrix(cur_col_points, rp_i, rn_i, r2_i, r3_i, SIE_quads, VIE_quads,...
                         res, 1, emu.k0, 1);
    
    rhs_E(:,i_sie) = rhs_E_i;
    rhs_H(:,i_sie) = rhs_H_i;
                                      
end


Ic = A \ rhs_E;
Ih = B \ rhs_H;

for i_sie = 1:N_sie

    count_vie = 1;
    j_vie = zeros(N_cells * PWX_sz,1);
    k_vie = zeros(N_cells * PWX_sz,1);
    cel_vox_idx = Cells(:, i_sie);    % returns array of global indices of cell voxels


    % loop over x,y,z directions
    for q = 1:dims.q
        % loop over basis component(s)
        for basis = 1:N_basis
            % loop over voxels in the cell
            for vox_num = 1:N_cells
                
                % get index for current voxel
                vox_idx = cel_vox_idx(vox_num);
                
                % form current indices of weights of expansion voxels
                j_vie(count_vie) = vox_idx + Next * (basis - 1) + Next * N_basis * (q - 1);
                k_vie(count_vie) = vox_num + N_cells * (basis-1) + N_cells * N_basis * (q - 1);
                
                count_vie = count_vie + 1;
                
            end
        end
    end
    
    Val{i_sie}   = Ic(k_vie, i_sie);
    I_val{i_sie} = j_vie;

end

% 
Val = vertcat(Val{:});
I_val = vertcat(I_val{:});

P = sparse(I_val, (1:N_sie) .* ones(N_cells * PWX_sz, 1), Val, Next_dofs, N_sie);


%% form mapping matrix from VIE to extended domain

% get expansion dofs indices
[~, I_Jb_idx, J_Jb_idx] = src_utils.get_exp_dofs(Cells, idxS_ext.', PWX_sz, Next);

V_Jb_idx = ones(size(I_Jb_idx));

S = sparse(I_Jb_idx, J_Jb_idx, V_Jb_idx, Next_dofs, N_vie_dofs);


% concatinate projector matrices
PS = [P,S];
