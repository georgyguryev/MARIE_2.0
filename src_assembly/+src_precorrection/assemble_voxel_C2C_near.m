function [Zc_fft] = assemble_voxel_C2C_near(mvp, projector, dims, freq, gpu_flag) 
% function compute_c2c_vox_representation(fL_op, COIL, RHBM, PROJ, freq)
% assembles coil-to-coil impedance matrix using voxelized representation

%% set parameters

% set electromagnetic parameters 
emu = src_utils.EM_utils(freq);

% get problem dimensions
N_sie  = dims.N_sie;
N_vie  = dims.N_scat * dims.ql;
N_ext  = dims.Nvox_ext;
N_near_1D = dims.N_near_1D;
N_near_3D = dims.N_near_3D;
N_exp_1D  = dims.N_exp_1D;
N_ext_dofs  = dims.Nvox_ext * dims.ql;
N_near_dofs = N_near_3D * dims.ql;

%% pre-allocate resulting matrix

Zc_fft = cell(N_sie,1);
I_fft  = cell(N_sie,1);
J_fft  = cell(N_sie,1);

Jcb = [speye(N_sie); sparse(N_vie, N_sie)];

% get coil 2 coil near lists
exp_cell_near = projector.cell_id;
c2c_near_list = projector.C2C_near_list;
near_vox_list = projector.near_vox_list;

% define 
exp_center      = floor(N_near_1D / 2) + mod(N_exp_1D, 2);
idx_exp_near_1D = src_projector.assign_expansion_cell(exp_center, exp_center, exp_center,...
                                                      N_exp_1D, N_near_1D,...
                                                      N_near_1D, N_near_1D); 
                                                  
% index exp near
idx_ql = src_scatterer.Scatterer_index(idx_exp_near_1D.', N_near_3D);
i_exp_near = repmat(idx_ql.index_ql(dims.ql, N_near_3D), 1, N_sie);
% i_exp_near = repmat(src_scatterer.Scatterer_index(idx_exp_near_1D.', N_near_3D).S_3d, 1, N_sie);

%% select appropriate L operator; allocate memory for buffer variable

if gpu_flag
    L_near    = @(Jin) mvp.L_near_gpu(Jin);
    Jc_near_i = zeros(N_near_dofs,1,'gpuArray');
else
    L_near    = @(Jin) mvp.L_near(Jin);
    Jc_near_i = zeros(N_near_dofs,1);
end

%% 

idx_sie    = 1:c2c_near_list.Count;
exp_voxels = exp_cell_near(:,idx_sie);
i_exp_ext  = src_scatterer.Scatterer_index(exp_voxels, N_ext).index_ql(dims.ql, N_ext);
j_exp_ext  = ones(size(i_exp_ext,1),1) * (1:N_sie);

% convert 2D subindices to 1D index
idx_exp_ext = sub2ind([N_ext_dofs N_sie], i_exp_ext(:), j_exp_ext(:));

% idx_near_1D = find(near_list);
Jc_ext  = mvp.Jcb2Jtot(Jcb);
Jc_near = sparse(i_exp_near, j_exp_ext, Jc_ext(idx_exp_ext), N_near_dofs, N_sie);



       
for i_sie = 1:c2c_near_list.Count
    
    near_list   = near_vox_list.ext(:,i_sie);
    idx_near_1D = near_vox_list.near{i_sie};
    c2c_near    = c2c_near_list(num2str(i_sie));
    
    near_list_dofs = src_scatterer.Scatterer_index(nonzeros(near_list), N_ext).index_ql(dims.ql, N_ext);
    idx_near_3D    = src_scatterer.Scatterer_index(idx_near_1D,N_near_3D).index_ql(dims.ql, N_near_3D);

    Jc_near_i(:) = full(Jc_near(:,i_sie));
    
    % apply Lop to projected unitary currents
    Ec_near = L_near(Jc_near_i);
    
    Ec_near = gather(Ec_near);
    
    % map near domain back to extended domain
    Ec_ext = sparse(near_list_dofs, ones(size(near_list_dofs)),...
                    Ec_near(idx_near_3D), N_ext_dofs, 1);
        
    % extract resulting fields on the coil from representation voxels
    Ecb = mvp.Jtot2Jcb(Ec_ext);
    
    Ecb = gather(Ecb);

    Zc_fft{i_sie} = Ecb(c2c_near);
    I_fft{i_sie}  = c2c_near';
    J_fft{i_sie}  = double(i_sie) * ones(size(c2c_near')); 
        
    clear Ec_near;

end

clear Jc_near_i;


% collapse all cell vectors to a single vector
Zc_fft = vertcat(Zc_fft{:});
I_fft  = vertcat(I_fft{:});
J_fft  = vertcat(J_fft{:});

% form a resulting sparse matrix for c2c precorrection
Zc_fft = sparse(I_fft, J_fft, Zc_fft, N_sie, N_sie);

Zc_fft = 1 / emu.ce * Zc_fft;
