function [Zc2b_Nop, Zc2b_Kop] = assemble_direct_B2C(coil, scatterer, projector, task_settings, dims, freq)
% function [Zbc_rwg] = assemble_direct_C2B(coil, scatterer, projector, task_settings, dims, freq)

%% setup parameters

% electromagnetic constants and frequency
emu = src_utils.EM_utils(freq);

% get mapping from scatterer index to scatterer
% degrees of freedom (non-air) voxels 
idxS_ext = scatterer.index_ext.S_1d;

% Get list of near scatterer voxels and 
C2B_near_list = projector.C2B_near_list;

% get size of near list
Nb     = dims.N_scat * dims.ql;
N_sie  = dims.N_sie;
N_scat = dims.N_scat;


% resolution 
res  = scatterer.dom_ext.res;  

% x,y,z coordinates of the domain
xd = scatterer.dom_ext.x_tensor;
yd = scatterer.dom_ext.y_tensor;
zd = scatterer.dom_ext.z_tensor;

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

%% compute near interactions between RWG and tissue voxels

% define cell arrays for Value, row and column indices respectively
Zbc_Nop_val = cell(C2B_near_list.Count,1);
Zbc_Kop_val = cell(C2B_near_list.Count,1);
I_idx   = cell(C2B_near_list.Count,1);
J_idx   = cell(C2B_near_list.Count,1);

for key = C2B_near_list.keys
    
     % current RWG basis number 
    i_sie = str2double(key{:});

    % get indices of scatterer voxels that are within RWG near domain
    near_scat_vox_ext = C2B_near_list(key{:});
    
    % define coordinated for scatterer voxels within RWG near domain
    Scoord = [xd(near_scat_vox_ext),...
        yd(near_scat_vox_ext),...
        zd(near_scat_vox_ext)];
    
    % find mapping from global scatterer indices to VIE dofs 1D
    [~, index_scat2dofs] = ismember(near_scat_vox_ext,idxS_ext);
    
    % generate dims.ql-dimensional VIE dofs indices
    idx_B_dofs = src_scatterer.Scatterer_index(index_scat2dofs, N_scat);
        
    % get coordinates of current rwg    
    rp_i = coil.rp(:,i_sie);
    rn_i = coil.rn(:,i_sie);
    r2_i = coil.r2(:,i_sie);
    r3_i = coil.r3(:,i_sie);
    
    % get fields produced by unitary excitation of current rwg triangle
    [Zbc_Nop_i, Zbc_Kop_i] = Assemble_rwg_coupling_matrix(Scoord.', rp_i, rn_i, r2_i, r3_i, SIE_quads, VIE_quads,...
                         res, dims.l, emu.k0, 1);
    
    Zbc_Nop_val{i_sie} = Zbc_Nop_i;
    Zbc_Kop_val{i_sie} = Zbc_Kop_i;
    I_idx{i_sie}   = idx_B_dofs.index_ql(dims.ql,dims.N_scat);
    J_idx{i_sie}   = i_sie * ones(size(Zbc_Nop_i));                                  
end



Zbc_Nop_val = vertcat(Zbc_Nop_val{:});
Zbc_Kop_val = vertcat(Zbc_Kop_val{:});
I_idx   = vertcat(I_idx{:});
J_idx   = vertcat(J_idx{:});

% form sparse coupling matrix
Zc2b_Nop = sparse(I_idx, J_idx, Zbc_Nop_val, Nb, N_sie);
Zc2b_Kop = sparse(I_idx, J_idx, Zbc_Kop_val, Nb, N_sie);
