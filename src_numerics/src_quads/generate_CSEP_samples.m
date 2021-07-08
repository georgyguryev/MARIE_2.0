function [r_src] = generate_CSEP_samples(basis, N_samples)
% function generates samples on a cylindrical equivalent surface 
% 
% Input parameters:
%
%           cyl_param - structure of parameters of an equivalent
%                       cylindrical surface 
%                       cyl_param.height - lenght of the cylinder
%                       cyl_param.rad    - radius of the cylinder
%
%           N_samples  - number of random samples, i.e. number of 
%                       elementary current sources (3D dipols) 
% Output parameters:
%           
%           r_sources - location of produced excitation sources 

%%

% get parameters
% define cylinder parameters
radius = basis.Cylinder_radius;
height = basis.Cylinder_height;

% get angle samples
phi_samp = 2 * pi * rand(N_samples,1);

% get hieght samples
height_samp = -height / 2  + height * rand(N_samples,1);

r_src = [radius * cos(phi_samp), radius * sin(phi_samp), height_samp];