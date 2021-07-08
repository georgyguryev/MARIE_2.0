function [Mo, Me, No, Ne] = calc_MNpot_bji(theta, phi, r, k, n)
% calculate the vector M N potential for sphere internal field while the Bessel
% function is the first kind Bessel function
% --input: theta, phi, r are the target position, they are scalars
%          k is the wavenumber, n is a vector of degrees
% --output: Mo, Me, No, Ne: 3-by-length(n) array
%          the 3 rows are for r, theta, phi directions
rho = k*r;

[~,~,bji, dbji, ~,~] = calc_sphbessel(n, 0, rho);

[Mo, Me, No, Ne] = calc_MNpot(theta, phi, rho, n, bji, dbji);
