function [RHBM] = properties_to_RHBM(r, rho, epsilon_r, sigma_e)
% function converts RHBM struct to arrays of parameters

if nargin == 4

    RHBM.r         = r;
    RHBM.rho       = rho;
    RHBM.epsilon_r = epsilon_r;
    RHBM.sigma_e   = sigma_e;
else
    RHBM = r;
end
