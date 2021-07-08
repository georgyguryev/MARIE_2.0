function [Npg,wt,Z,z1,z2,z3] = gauss_2d_triangles(Np)

format long

% get 1D gauss quadrature points
[w,z] = gauss_1d (Np);

% compute auxiliary variables
W  = w * w';
x  = kron(eye(Np), ones(Np,1)) * z; 
y  = repmat(z,Np,1);
zi = (1 - y) / 8;

% quadrature weights
wt = (W(:).*zi)'; 

% quadrature points
z2     = (1 + y) / 2;
z3     = (1 - z2) .* (1 + x) / 2;
z1     = 1 - z2 - z3;

% dimensions of 2D Gauss quadratures 
Npg = Np * Np;

% concatenate quadrature points to [Np x 3] matrix Z
Z=[z1 z2 z3];
