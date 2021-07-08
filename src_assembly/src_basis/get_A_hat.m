function [A_hat_out] = get_A_hat(rows, cols, task_settings, scatterer, coil, dims, freq, op_type)
%  function get_A_hat() assembles columns of coupling operator; 
%  This function is written to be used in cross_2d() to addaptively 
%  construct a cross-copressed representation of the dense coupling operator
%   Input parameters: 
%       cols          - RWG-related columns to be sampled 
%       task_settings - contains all coupling-related settings
%       dims          - an object that contains all problem-specific dimensions
%       scatterer     - contains a scatterer object (mat. properties,
%                       masks, etc.
%       coil          - contains all coil related data
%       freq          - current working frequency
%       op_type       - 'N' for N and 'K' for K operators; 
%                       assembles both otherwise

A_hat_N = [];
A_hat_K = [];

[vox,l,q] = ind2sub([dims.N_scat, dims.l, dims.q],  rows);
[I_vox, ~, ic] = unique(vox);

i_mask =  sub2ind([size(I_vox,1), dims.l, dims.q], ic,l,q);
I = speye(dims.ql * size(I_vox,1));
P_mask = I(i_mask,:);


switch op_type
    case 'N'                                                
        [A_hat_N, ~] = src_coupling.assemble_coupling_matrices(task_settings, scatterer, coil, dims, freq, cols, I_vox, 0);
        A_hat_N = P_mask * A_hat_N;
    case 'K'
        [~, A_hat_K] = src_coupling.assemble_coupling_matrices(task_settings, scatterer, coil, dims, freq, cols, I_vox, 0);
        A_hat_K = P_mask * A_hat_K;

    otherwise 
        error('Cross: Error! get_A_hat() supports only  "N" and "K" operators');
end

A_hat_out = [A_hat_N, A_hat_K];
% EOF
end