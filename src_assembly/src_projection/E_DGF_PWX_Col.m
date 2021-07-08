function    [Nop, Kop] =  E_DGF_PWX_Col(r_center, res, r_col, GL_order, freq, dims) 
% function  [E] =  E_DGF_PWX_Col(r_center, res, r_col, GL_order, freq, dims) 
% this function computes three components of the electric fields 
% created by expansion voxels at collocation locations r_col
% -----------------------------------------------------------------------
%                      INPUT:
% -----------------------------------------------------------------------
%       r_center - coordinates of centers of expansion voxels
%       res      - resolution
%       r_col    - coordinates of centers of collocation points
%       GL_order - order of 1D Gauss quadrature (VIE)
%       freq     - frequency 
%       dims     - structure that keeps dimensions of all objects    

%% Updated 04/19/2021 by Georgy Guryev, RLE

%% get size of output vector
N_basis = dims.l;

% get EM constants
emu = src_utils.EM_utils(freq);

% get wave number
ko     = emu.k0;
%% get 1D VIE points 

[wt_vie,z] = gauss_1d(GL_order);
vie_quad   = [size(wt_vie,1); wt_vie; z];

%% this part is developed specificaly for a single PWC sopport

% construct operator that describes coupling between expansion voxels and
% collocation points
[Nop,Kop] = Assemble_p2v_coupling(r_center, r_col.', vie_quad, res, N_basis, ko);

% reshape matrix to the appropriate form
N_dofs = size(Nop,1);
Nop      = permute(reshape(Nop.', 3, [], N_dofs), [2 1 3]);
Kop      = permute(reshape(Kop.', 3, [], N_dofs), [2 1 3]);
Nop      = res^3 * reshape(Nop, [], N_dofs);
Kop      = res^3 * reshape(Kop, [], N_dofs);
