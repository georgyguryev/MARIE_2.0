function Vout = coupling_c2b_basis_mvp(Jin, U, alpha, Gram, S_q2ql)
% function coupling_c2b_basis_mvp() computes fields 
% incident on the scatterer
% INPUT: 
%       - Jin  - volumetric tissue currents 
%       - U  - left sub-space of the coupling operator
%       - alpha - coil-specific basis weights
%       - Gram  - basis Gramian
%       - S_q2ql - PWC to PWL mask 
% OUTPUT:
%        - Vout - electric fields, incident on the tissue
% Author: Georgy Guryev, Cambridge, MA, 2019   


% if nargin == 5  && size(U,1) ~= size(S_q2ql,1)
%     Vout = Gram * S_q2ql * (U * (alpha * Jin));
% else 
%     Vout = Gram * (U * (alpha * Jin));
% end

if nargin == 5  && size(U,1) ~= size(S_q2ql,1)
    Vout = Gram * S_q2ql * (U * (alpha * Jin));
else 
    Vout = Gram * (U * (alpha * Jin));
end