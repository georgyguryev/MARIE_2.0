function [U_out] = get_U(cols, task_settings, scatterer, coil, dims, freq, op_type)
%  function get_U() assembles columns of coupling operator; 
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

U_N = [];
U_K = [];

switch op_type
    case 'N'
        [U_N, ~] = src_coupling.assemble_coupling_matrices(task_settings, scatterer, coil, dims, freq, cols, [], 0);
    case 'K'
        [~, U_K] = src_coupling.assemble_coupling_matrices(task_settings, scatterer, coil, dims, freq, cols, [], 0);
    otherwise
        [U_N, U_K] = src_coupling.assemble_coupling_matrices(task_settings, scatterer, coil,  dims, freq, cols, [], 0);
end

U_out = [U_N, U_K];

% EOF
end

