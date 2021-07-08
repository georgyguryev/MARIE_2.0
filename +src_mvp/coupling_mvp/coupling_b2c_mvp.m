function Vout = coupling_b2c_mvp(Jin, mvp, precorrection, dims, freq)
% function coupling_b2c_mvp() computes fields 
% incident on the scatterer via pFFT
% INPUT: 
%       - Jin  - volumetric tissue currents 
%       - mvp  - matrix-vector products
%       - precorrection - coupling precorrection
%       - dims  - object with problem dimensions
%       - freq - working frequency 
% OUTPUT:
%        - Vout - electric fields, incident on the coil 
% Author: Georgy Guryev, Cambridge, MA, 2019    

% get elecromagnetic constants
emu = src_utils.EM_utils(freq);

%% for Jcb vector from Jin == Jb

Jcb = [zeros(dims.N_sie, 1); Jin];

%% compute total current and map volumetric current

% map Jc, Jb to J_tot 
J_tot = mvp.Jcb2Jtot(Jcb);

%% compute L operator for L(Jb_ext) and for L(Jtot)

% Compute fields, produced by coil surface currents within extended domain
Vb_vox = mvp.L_ext(J_tot);

% get fields within a scatterer domain 
Vcb = mvp.Jtot2Jcb(Vb_vox(:));

% scale field intensities with electromagnetic constants
Vout = 1 / emu.ce * Vcb(1:dims.N_sie) + precorrection.Z_B2C.' * Jin;