function [V_out] = get_V_bb(rows, task_settings, scatterer, coil, dims, freq, op_type)
%  function get_V_bb() assembles rows of coupling operator; 
%  This function is written to be used in cross_2d() to addaptively 
%  construct a cross-copressed representation of the dense coupling operator
%   Input parameters: 
%       rows          - voxel-related rows to be sampled 
%       task_settings - contains all coupling-related settings
%       dims          - an object that contains all problem-specific dimensions
%       scatterer     - contains a scatterer object (mat. properties,
%                       masks, etc.
%       coil          - contains all coil related data
%       freq          - current working frequency
%       op_type       - 'N' for N and 'K' for K operators; 
%                       assembles both otherwise

V_N = [];
V_K = [];

% find voxels that should be sampled
[vox,l,q] = ind2sub([dims.Nvox_vie, dims.l, dims.q],  rows);
[I_vox, ~, ic] = unique(vox);

% filter out redundant components 
i_mask =  sub2ind([size(I_vox,1), dims.l, dims.q], ic,l,q);
I = speye(dims.ql * size(I_vox,1));
P_mask = I(i_mask,:);


switch op_type
    case 'N'                                                
        [V_N, ~] = src_coupling.assemble_coupling_matrices(task_settings, scatterer, coil, dims, freq, [], I_vox, 1);
        V_N = P_mask * V_N;
    case 'K'
        [~, V_K] = src_coupling.assemble_coupling_matrices(task_settings, scatterer, coil, dims, freq, [], I_vox, 1);
        V_K = P_mask * V_K;
    otherwise
        [V_N, V_K] = src_coupling.assemble_coupling_matrices(task_settings, scatterer, coil,  dims, freq, [], I_vox, 1);
        V_N = P_mask * V_N;
        V_K = P_mask * V_K;
end

V_out = [V_N; V_K];

% EOF
end
