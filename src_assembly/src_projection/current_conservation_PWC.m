function [I] = current_conservation_PWC(r_center, res, GL_order)

[m,n] = size(r_center);

% Allocte memory 
I = zeros(m,m*n);

% define indices of non-zero components
idx_x = 1:3:m*n-2;
idx_y = 2:3:m*n-1;
idx_z = 3:3:m*n;

% define values of non-zero current components
I(1,idx_x) = res^3;
I(2,idx_y) = res^3;
I(3,idx_z) = res^3;
