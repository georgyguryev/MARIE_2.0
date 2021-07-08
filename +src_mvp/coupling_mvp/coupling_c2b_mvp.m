function Vout = coupling_c2b_mvp(Jin, mvp, precorrection, dims, freq)
% function coupling_c2b_mvp() computes a matrix-vector product between
% implicit couping operator and input coil currents; the function returns
% the resulting incident fields within body domain
% INPUT:
%       - Jin - surface coil currents
%       - mvp - matrix-vector products
%       - precorrection - coupling precorrection
%       - dims - object with problem dimensions
%       - freq - working frequency
% OUTPUT:
%       - Vout - fields, incident on the tissue
% Author: Georgy Guryev, Cambridge, MA, 2019   

% get elecromagnetic constants
emu = src_utils.EM_utils(freq);

%% for Jcb vector from Jin == Jc

Jcb = [Jin; zeros(dims.N_scat * dims.ql, 1)];

%% compute total current and map volumetric current

% map Jc, Jb to J_tot 
J_tot = mvp.Jcb2Jtot(Jcb);

%% compute L operator for L(Jb_ext) and for L(Jtot)

% Compute fields, produced by coil surface currents within extended domain
Vc_vox = mvp.L_ext(J_tot);

% get fields within a scatterer domain 
Vcb = mvp.Jtot2Jcb(Vc_vox(:));

% scale field intensities with electromagnetic constants
Vout = 1 / emu.ce * Vcb(dims.N_sie + 1:end) + precorrection.Z_C2B.Nop * Jin;