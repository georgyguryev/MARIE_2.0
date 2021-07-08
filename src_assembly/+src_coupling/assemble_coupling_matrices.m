function [Zbc_N, Zbc_K] = assemble_coupling_matrices(task_settings, scatterer, coil, dims, freq, mask, samp_vox, flag_bb)
% function assemble_coupling_matrices(task_settings, scatterer, coil, dims, freq, mask, samp_vox, flag_bb)
% simultaneously assembles electric and magnetic coupling operators (or
% their subsets)
% INPUT:
%        - task_settings - simulation setup from .ini file
%        - scatterer     - scatterer object
%        - coil          - coil object
%        - dims          - object stores all dimensions
%        - mask          - 
%        - freq          - simulation frequency
%        - samp_vox      - sampled tissue voxels
%        - flag_bb       - flag indicates if the coupling operator should
%        be assembled for the tissue or the bounding box, containing one
% OUTPUT:
%        - Zbc_N - coupling operator for the electric field
%        - Zbc_K - coupling operator for the magnetic field
% Author: Georgy Guryev, Cambridge, MA, 2019    


l   = dims.l;
emu = src_utils.EM_utils(freq);
k0  = emu.k0;


%
Quad_order_vie = task_settings.vsie.Np_quad_coup_vie;
Quad_order_sie = task_settings.vsie.Np_quad_coup_sie;

[wp_vie, z_vie] = gauss_1d(Quad_order_vie);
[Np_sie, z1_sie, z2_sie, z3_sie, wp_sie] = dunavant_rule(Quad_order_sie);

VIE_quads = [Quad_order_vie; wp_vie; z_vie];
SIE_quads = [Np_sie; wp_sie; z1_sie; z2_sie; z3_sie];

% x,y,z coordinates of the domain
xd = scatterer.dom_vie.x_tensor;
yd = scatterer.dom_vie.y_tensor;
zd = scatterer.dom_vie.z_tensor;
idxS = scatterer.index_vie.S_1d;
res =  scatterer.dom_vie.res;
vol = res^3;


if flag_bb == 1
    idxS = (1:dims.Nvox_vie).';
end

if ~isempty(samp_vox)
    if flag_bb == 1
        idxS = samp_vox;
    else
        idxS = idxS(samp_vox);
    end
end

Scoord = [xd(idxS) yd(idxS) zd(idxS)].';

if nargin < 6 || isempty(mask)
    mask = 1:dims.N_sie;
end

rp = coil.rp(:,mask);
rn = coil.rn(:,mask);
r2 = coil.r2(:,mask);
r3 = coil.r3(:,mask);

% assemble matrix for specified triangles
[Zbc_N, Zbc_K] = Assemble_rwg_coupling_matrix(Scoord, rp, rn, r2, r3, SIE_quads, VIE_quads,...
    res, l, k0, 1);

% scale coupling matrices with voxel volume
Zbc_N = Zbc_N * vol;
Zbc_K = Zbc_K * vol;