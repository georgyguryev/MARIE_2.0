function [Zc_fft] = assemble_voxel_C2C(mvp, dims, freq, gpu_flag) 
% function compute_c2c_vox_representation(fL_op, COIL, RHBM, PROJ, freq)
% assembles coil-to-coil impedance matrix using voxelized representation

%% set parameters

% set electromagnetic parameters 
emu = src_utils.EM_utils(freq);

% get problem dimensions
N_sie = dims.N_sie;
N_vie = dims.N_scat * dims.ql;

%% pre-allocate resulting matrix
Zc_fft = zeros(N_sie);


Jcb = [eye(N_sie); zeros(N_vie, N_sie)];

% select appropriate L operator; allocate memory for buffer variable
if gpu_flag
    
    Jcb2Jtot    = @(Jcb)  mvp.Jcb2Jtot_gpu(Jcb);
    Jtot2Jcb    = @(Jtot) mvp.Jtot2Jcb_gpu(Jtot);
    L_ext       = @(Jin)  mvp.L_ext_gpu(Jin);
    
    Jcb = gpuArray(Jcb);
else
    Jcb2Jtot    = @(Jcb)  mvp.Jcb2Jtot(Jcb);
    Jtot2Jcb    = @(Jtot) mvp.Jtot2Jcb(Jtot);
    L_ext       = @(Jin)  mvp.L_ext(Jin);
end


%% 
tic
for i_sie = 1:N_sie
    
    Jcb_ii = Jcb(:,i_sie);
%     tic
    Jc_ext = Jcb2Jtot(Jcb_ii); 
       
    % apply Lop to projected unitary currents
    Ec_ext_loc = L_ext(Jc_ext);
        
    % extract resulting fields on the coil from representation voxels
    Ecb = 1 /emu.ce * Jtot2Jcb(Ec_ext_loc(:));
    
    Zc_fft(:,i_sie) = gather(Ecb(1:N_sie));
    
    clear Ecb Ec_ext_loc Jcb_ii Jc_ext;
%     toc
end
toc
