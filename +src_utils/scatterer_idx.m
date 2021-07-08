function idx = scatterer_idx(RHBM, freq)
% fucntion returns indices of scatterer voxels

if nargin < 2
    freq = 298e6;
end

epsilon_r = RHBM.epsilon_r;
sigma_e   = RHBM.sigma_e;

emu = src_utils.EM_utils(freq);

e_r = epsilon_r + 1i * sigma_e / emu.ce;

% find non-air voxels (check if e_r is below the threshold
idx = find(abs(e_r(:) - 1) > 1e-12);