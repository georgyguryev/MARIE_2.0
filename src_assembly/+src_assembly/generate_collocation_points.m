function [col_points] = generate_collocation_points(res, col_distance, N_col_pts)
% function generate_collocation_points() generates collocation points on an
% sphere of radius col_distance using LebedevSphere
% INPUT: 
%        - res - residual
%        - col_distance  - distance of the collocation sphere
%        - N_coil_pts    - degree of collocation points
% OUTPUT:
%        - col_points - coordinates of collocation points
% Author: Georgy Guryev, Cambridge, MA, 2019

radius = col_distance .* res;

% generate Lebedev Quadratures
leb = getLebedevSphere(N_col_pts);

col_points = radius .* [leb.x.'; leb.y.'; leb.z.'];


