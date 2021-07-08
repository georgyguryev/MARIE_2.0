function [Hout] = Jc2H_coupling_basis(Jin, U_K, alpha_K, S_q2ql)
% function Vout = coupling_c2b_basis_mvp(Jin, U, alpha, Gram, S_q2ql)
% function coupling_c2b_basis_mvp() computes fields 
% incident on the scatterer (for now assumes that BB is scatterer)

% if nargin == 5  && size(U,1) ~= size(S_q2ql,1)
%     Vout = Gram * S_q2ql * (U * (alpha * Jin));
% else 
%     Vout = Gram * (U * (alpha * Jin));
% end

if nargin == 5  && size(U_K,1) ~= size(S_q2ql,1)
    Hout = S_q2ql * (U_K * (alpha_K * Jin));
else 
    Hout = U_K * (alpha_K * Jin);
end