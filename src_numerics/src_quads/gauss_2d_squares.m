function [Npg,w,u,v] = gauss_2d_squares(Np)

% get 1D gauss quadrature points
[w_1d, x_1d] = gauss_1d(Np);

% define number of 2D quadrature points
Npg = Np  * Np;

% define new weights and quadrature points for square domain ||x||_Inf = 1
w = reshape(w_1d * w_1d',[],1) / 4; 
u = (kron(eye(Np), ones(Np,1)) * x_1d + 1) / 2; 
v = (repmat(x_1d, Np, 1) + 1) / 2;
