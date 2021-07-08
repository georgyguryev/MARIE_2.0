function [J_hat, rank_hat] = sample_cols(V, J, rank, N, tol)

% leb_degrees = [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, ...
%              350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, ...
%              3470, 3890, 4334, 4802, 5294, 5810];

if isempty(V)
    J_hat = sort(randi(N, rank, 1));
else
    J_hat = sort(maxvol(V, tol, rank^2));
    J_hat = setdiff(J_hat, J);
%     [~,i_ld] = min(abs(leb_degree - N));
%     degree = leb_degrees(i_ld);
%     leb_quad = getLebedevSphere(degree);
%     
%     rxy  = sqrt(leb_quad.x.^2 + leb_quad.y.^2);
%     
%     theta = acos(leb_quad.z);
%     phi   = acos(leb_quad.x ./ rxy);
%     JJ_hat = 
end

rank_hat = size(J_hat,1);