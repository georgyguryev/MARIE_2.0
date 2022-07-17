function Vout = coupling_UWU_basis_mvp(Jin, U, W, M_q2ql)
% function implements modified basis-based mvp 

if nargin == 4  && size(U,1) ~= size(M_q2ql,1)
    Vout =  M_q2ql * (U * (W  * (M_q2ql.' * Jin)));
else 
    Vout = U * (W  * Jin);
end

% Vout = gather(Vout);

