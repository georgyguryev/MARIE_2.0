function [E_bound] = mvp_P(Jexp, fN_mvp_near, projector, dims)
% function computes fields, produced by fictitious volumetric currents at
% the boundary of 'near zone'; The fields at the boundary are used to
% assemble the Projection matrix   

% get problem dimensions 
N_near_1d   = dims.N_near_1D;
N_exp_1d    = dims.N_exp_1D;
frame_width = projector.near_boundary_width;

% for collocation at the boundary voxels
in_near_st  = frame_width + 1;
in_near_end = N_near_1d - frame_width;

% get number of basis functions
PWX_sz   = dims.ql;

% allocate Jin = [N_near x N_near x N_near] 
Jin = zeros(N_near_1d, N_near_1d, N_near_1d, PWX_sz);

% find shift from left bottom corner of near domain to exp. domain
idx_exp_shift = round((N_near_1d - N_exp_1d) / 2);

% define range for expansion domain voxels (due to simmatry idx_y and z are the same)
idx_exp_x = idx_exp_shift + (1:N_exp_1d);

% % reshape vector Jexp to tensor
Jexp = squeeze(reshape(Jexp, dims.exp));

% copy nonzero components to Jin
Jin(idx_exp_x,idx_exp_x,idx_exp_x,:) = Jexp;

% apply N operator 
E_near = fN_mvp_near(Jin(:));

% reshape and set to zero all non-boundary field components
E_near = reshape(E_near, N_near_1d, N_near_1d, N_near_1d, PWX_sz);
E_near(in_near_st:in_near_end, in_near_st:in_near_end, in_near_st:in_near_end, :) = 0;

%% get fields at the boundary; return results

E_bound = E_near(E_near~=0);

% E_bound = squeeze(E_near(N_near_1d,N_near_1d,N_near_1d,:));

