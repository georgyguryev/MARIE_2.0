function [E] = fill_E_DGF_PWX(E_cur, rows, cols, i_x, i_y, i_z, j_x, j_y, j_z)
% ----------------------------------------------------------------------- %
% function   [E] = fill_E_DGF_PWX(E,i_str, i_end, j_x,j_y,j_z)
%   fills entries of matrix E (electric field components at col. points)
%   The function was written to enable parfor in function E_DGF_PWX_col();
%
% ----------------------------------------------------------------------- %
%   INPUT:
%   
%           E_cur - current field components to be inserted
%           rows, cols - dimentions of matrix to be filled
%           i_str, i_end - range of rows to be filled
%           j_x,j_y,j_z - indicies of voxels to be filled
% ----------------------------------------------------------------------- %

%% 

E = zeros(rows,cols);


% three field components produced by tissue current along x-axis
E(i_x,j_x) = E_cur(1,1);
E(i_y,j_x) = E_cur(2,1);
E(i_z,j_x) = E_cur(3,1);

% three field components produced by tissue current along y-axis
E(i_x,j_y) = E_cur(1,2);
E(i_y,j_y) = E_cur(2,2);
E(i_z,j_y) = E_cur(3,2);


% three field components produced by tissue current along y-axis
E(i_x,j_z) = E_cur(1,3);
E(i_y,j_z) = E_cur(2,3);
E(i_z,j_z) = E_cur(3,3);
