function [Zbc_Nop, Zbc_Kop] = assemble_voxel_B2C(scatterer,projector,mvp, dims, gpu_flag)
% -------------------------------------------------------------------------
% function [Zbc_fft] = compute_coup_near_FFT(RHBM,SCOIL,PROJ,freq,gpu_flag)
%
% this fuction computes interactions between representation voxels and body
% voxels that lie within the near zone of a given RWG. This function has to
% be performed for both implicit coupling and full pFFT implementations
% -------------------------------------------------------------------------
%
% Input parameters:
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
%
% Output parameters:
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------

%% setup parameters


idxS_ext = scatterer.index_ext.S_1d;
% idx_exp  = projector.index_exp2near.S_3d;
idx_exp  = projector.index_exp2near.index_ql(dims.ql, dims.N_near_3D);


% list of near interactions
pwx_near_list = projector.near_vox_list.ext;
C2B_near_list = projector.C2B_near_list;

% get size of near list
N_sie       = dims.N_sie;
N_scat      = dims.N_scat;
N_exp       = dims.N_exp_3D;
N_exp_ql    = dims.ql * N_exp;
N_near      = dims.N_near_3D;
N_near_ql   = dims.ql * N_near;
N_ext       = dims.Nvox_ext;
Nb          = dims.N_scat * dims.ql;

% exp2ext index
exp2ext_dofs_idx = repmat(projector.cell_id, dims.ql, 1) + ...
                   kron(N_ext * [0:dims.ql-1].', ones(N_exp, N_sie));
                   
% select coil-related part of the projection matrix 
P = projector.PS(:,1:N_sie);

% select appropriate L operator; allocate memory for buffer variable
if gpu_flag
    L_near       = @(Jin) mvp.L_near_gpu(Jin);
    Lnear_exp = gpuArray(zeros(N_near_ql, N_exp_ql));
    K_near    = @(Jin) mvp.K_near_gpu(Jin);
    Knear_exp = gpuArray(zeros(N_near_ql, N_exp_ql));
else
    L_near = @(Jin) mvp.L_near(Jin);
    Lnear_exp = zeros(N_near_ql, N_exp_ql);
    K_near    = @(Jin) mvp.K_near(Jin);
    Knear_exp = zeros(N_near_ql, N_exp_ql);
end



%% compute near interactions

% get triangles that have scatterer voxels in near zone
keys = C2B_near_list.keys;

Zbc_Nop_val = cell(length(keys),1);
Zbc_Kop_val = cell(length(keys),1);
I_idx   = Zbc_Nop_val;
J_idx   = I_idx;

P_hat = zeros(dims.ql * dims.N_exp_3D, length(keys));

for i = 1:length(keys)
    
    % current RWG basis number
    i_sie = str2double(keys{i});
    
    % get current expansion list
    exp_idx = exp2ext_dofs_idx(:,i_sie);
    P_hat(:,i) = P(exp_idx, i_sie);
end

I   = eye(N_exp_ql);
Jin = zeros(N_near_ql, N_exp_ql);
Jin(idx_exp,1:N_exp_ql) = I(:,1:N_exp_ql);

if gpu_flag
    Jin = gpuArray(Jin);
end


for j = 1:N_exp_ql
    Lnear_exp(:,j) = reshape(L_near(Jin(:,j)),[],1);
    Knear_exp(:,j) = reshape(K_near(Jin(:,j)),[],1);
end

Lnear_exp = gather(Lnear_exp);
Knear_exp = gather(Knear_exp);

E_near = Lnear_exp * P_hat;
K_near = Knear_exp * P_hat;

for i = 1:length(keys)
    
    % current RWG basis number 
    i_sie = str2double(keys{i});    
    
    % get indices of scatterer voxels that are within RWG near domain
    near_scat_vox_ext = C2B_near_list(keys{i});
    
    % get global indices of near voxels for current RWG
    near_vox_ext      = pwx_near_list(:,i_sie);
    
    % find local indices of scatterer voxels in 
    near_scat_vox_local  = find(ismember(near_vox_ext, near_scat_vox_ext)); 
    
    % find mapping from global scatterer indices to VIE dofs
    index_near_scat2dofs = find(ismember(idxS_ext, near_scat_vox_ext));
    
    % find local indices of scatterer degrees of freedom
    index_near_scat_loc  = src_scatterer.Scatterer_index(near_scat_vox_local, N_near);
    
    % find global indices of scatterer degrees of freedom
    index_near_scat_dofs = src_scatterer.Scatterer_index(index_near_scat2dofs, N_scat);
    

    % get fields at scatterer voxels and store them in a vector
    Zbc_Nop_val{i} = E_near(index_near_scat_loc.index_ql(dims.ql,N_near),i);
    Zbc_Kop_val{i} = K_near(index_near_scat_loc.index_ql(dims.ql,N_near),i);
    I_idx{i}   = index_near_scat_dofs.index_ql(dims.ql,dims.N_scat);
    J_idx{i}   = i_sie * ones(size(index_near_scat_loc.index_ql(dims.ql,dims.N_scat)));
    
end


% construct a sparse matrix from index vectors and a value vector
Zbc_Nop_val = vertcat(Zbc_Nop_val{:});
Zbc_Kop_val = vertcat(Zbc_Kop_val{:});

I_idx   = vertcat(I_idx{:});
J_idx   = vertcat(J_idx{:});

% form sparse coupling matrices 
Zbc_Nop = sparse(I_idx, J_idx, Zbc_Nop_val, Nb, N_sie);
Zbc_Kop = sparse(I_idx, J_idx, Zbc_Kop_val, Nb, N_sie);






