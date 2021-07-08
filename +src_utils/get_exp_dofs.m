function [dofs_idx, I_Jb_idx, J_Jb_idx] = get_exp_dofs(Cells,idxS_ext, PWX_sz, N_ext) 
% function get_exp_dofs(Cells, N_sie, N_exp, PWX_sz) maps global index
% mapping (extended domain) to dofs index space, i.e. from [L x M x N] to 
% [L x M x N x [l x q]]

% N_exp_dofs = N_exp * PWX_sz;
N_ext_dofs = N_ext * PWX_sz;

% % allocate memory for expension domain dofs index
% dofs_idx = zeros(N_exp_dofs, N_sie);

    
% get unshifted idx of expansion dofs
usft_dofs_idx    = repmat(Cells, PWX_sz, 1);
usft_dofs_Jb_idx = repmat(idxS_ext,1, PWX_sz);

% define shift scaling  
shift_factor    = 0:N_ext:N_ext_dofs - N_ext;

% form shift due to l and q components
shift_exp = kron(shift_factor.',ones(size(Cells)));
shift_Jb  = kron(shift_factor,ones(size(idxS_ext)));

% get idx of expansion dofs
dofs_idx = usft_dofs_idx + shift_exp;
I_Jb_idx = usft_dofs_Jb_idx + shift_Jb;
J_Jb_idx = 1:size(I_Jb_idx,2);
