function [r] = generate_extended_domain(Xext, res)
% function [r] = generate_extended_domain(X_ext, res)
% generates extended uniform grid for SCOIL+RHBM domain (pFFT)
% ------------------------------------------------------------------------
% INPUT :
% ------------------------------------------------------------------------
%     Xest - 3x2 matrix of domain limits (gives centers of boundary voxels)
%            Xext = [ xmin xmax; ymin ymax; zmin zmax]
%     res  - voxel resolution / length of voxel edge  
% ------------------------------------------------------------------------
% OUTPUT :
% ------------------------------------------------------------------------
%     r  - 4D  structure that stores grid coordinates
% ------------------------------------------------------------------------


% -------------------------------------------------------------------------
% generate coordinate arrays
% -------------------------------------------------------------------------

% obtain length of each size
% Dx = Xext(1,2) - Xext(1,1);
% Dy = Xext(2,2) - Xext(2,1);
% Dz = Xext(3,2) - Xext(3,1);

% obtain minimum number of cells in each direction
% Lext = round(Dx/res) + 1;
% Mext = round(Dy/res) + 1;
% Next = round(Dz/res) + 1;
% Lext = ceil(Dx/res);
% Mext = ceil(Dy/res);
% Next = ceil(Dz/res);


% fprintf('[Lext x Next x Mext] = %d %d %d \n', Lext, Next, Mext);


% generate the arrays
% x = linspace(Xext(1,1),Xext(1,2), Lext + 1);
% y = linspace(Xext(2,1),Xext(2,2), Mext + 1);
% z = linspace(Xext(3,1),Xext(3,2), Next + 1);

% x = Xext(1,1):res:Xext(1,1) + Lext * res;
% y = Xext(2,1):res:Xext(2,1) + Mext * res;
% z = Xext(3,1):res:Xext(3,1) + Next * res;

x = Xext(1,1):res:Xext(1,2);
y = Xext(2,1):res:Xext(2,2);
z = Xext(3,1):res:Xext(3,2);

% -------------------------------------------------------------------------
% Generate grid
% -------------------------------------------------------------------------

r = src_geo.grid3d(x,y,z);