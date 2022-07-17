function [Hout] = Jb2H_coupling_basis(Jin, Ic0, U_K, V_K, b_K, b_N, Y_LC, rhs_cp, S_q2ql)
% function Vout = coupling_c2b_basis_mvp(Jin, U, alpha, Gram, S_q2ql)
% function coupling_c2b_basis_mvp() computes fields 
% incident on the scatterer (for now assumes that BB is scatterer)

% if nargin == 5  && size(U,1) ~= size(S_q2ql,1)
%     Vout = Gram * S_q2ql * (U * (alpha * Jin));
% else 
%     Vout = Gram * (U * (alpha * Jin));
% end

Np = size(Y_LC, 2);
Ib = eye(Np);

if nargin == 9  && size(U_K,1) ~= size(S_q2ql,1)
    Hout = S_q2ql * (b_K * rhs_cp - U_K * (V_K * Jin));
else 
    Hout = b_K * (Ib - Y_LC \ Ic0) * rhs_cp - U_K * (V_K * Jin) - (b_K / Y_LC) * (b_N.' * Jin);
end