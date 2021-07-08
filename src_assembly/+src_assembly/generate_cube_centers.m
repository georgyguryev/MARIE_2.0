function [r_center] = generate_cube_centers(r_c, res, L, M, N)
% function [r_center] = generate_cube_centers(r_c, res, L, M, N)
% generates coordinates of centers of expansion voxels
% INPUT:
%        - r_c - coordinates of center of expansion domain
%        - res - resolution
%        - L, M, N - number of voxels along x,y and z axis
% OUTPUT:
%        - r_center - coordinates of centers of expansion voxels
% Author: Georgy Guryev, Cambridge, MA, 2019

if(nargin < 5)
    M = L;
    N = L;
end

% number of expansion voxels along each axis respectively
N_disp_x = floor(L / 2);
N_disp_y = floor(M / 2);
N_disp_z = floor(N / 2);

x = res * linspace(0,L - 1,L) - res * N_disp_x;
y = res * linspace(0,M - 1,M) - res * N_disp_y;
z = res * linspace(0,N - 1,N) - res * N_disp_z;

r_o = src_geo.grid3d(x,y,z);

% 
X_o = r_o(:,:,:,1) + r_c(1);
Y_o = r_o(:,:,:,2) + r_c(2);
Z_o = r_o(:,:,:,3) + r_c(3);

% allocate memory for tensor
r_center = zeros(size(r_o));

% copy voxel centers back to tensor
r_center(:,:,:,1) = X_o;
r_center(:,:,:,2) = Y_o;
r_center(:,:,:,3) = Z_o;

r_center = round(reshape(r_center,[L * M * N, 3]),4);
 


