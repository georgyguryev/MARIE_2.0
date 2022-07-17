function Vout = coupling_b2c_basis_mvp(Jin, U, alpha, Gram, S_q2ql)
% function coupling_b2c_basis_mvp() computes fields 
% incident on the scatterer 
% INPUT: 
%       - Jin  - volumetric tissue currents 
%       - U  - left sub-space of the coupling operator
%       - alpha - coil-specific basis weights
%       - Gram  - basis Gramian
%       - S_q2ql - PWC to PWL mask 
% OUTPUT:
%        - Vout - electric fields, incident on the coil 
% Author: Georgy Guryev, Cambridge, MA, 2019    

if nargin == 5  && size(U,1) ~= size(S_q2ql,1) && ~isempty(S_q2ql)
    Vout = alpha.' * (U.' * ( S_q2ql.' * Jin));
else
    Vout = alpha.' * (U.' * Jin);
end