function [r, rho, epsilon_r, sigma_e] = RHBM_to_properties(RHBM)
% function converts RHBM struct to arrays of parameters


r         = RHBM.r;
rho       = RHBM.rho;
epsilon_r = RHBM.epsilon_r;
sigma_e   = RHBM.sigma_e;
