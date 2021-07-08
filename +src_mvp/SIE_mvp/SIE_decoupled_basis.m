function Vout = SIE_decoupled_basis(Jb, feed_port, U, X, Jc_ini, M_q2ql)
% function [Vout] = SIE_decoupled_basis(Jb, feed_port, U, X, Jc_ini, M_q2ql)
% computes coil surface currents from inital guess (SIE solution) and
% contribution of tissue currents (basis)
% ------------------------------------------------------------------------
% Input:
%        - Jb - tissue polarization currents (vector)
%        - feed_port - index of the feed/excitation port
%        - U,X - object that contains all tissue specific information
%        - Jc_ini - initial solution for surface currents (self-consistant)
%        - M_q2ql - PWC/PWL masks
% Output:
%        - Vout - the resulting solution for surface currents 
% ------------------------------------------------------------------------ 
% Author: Georgy Guryev, Cambridge, MA, 2019    

if nargin == 6  && size(U,1) ~= size(M_q2ql,1)
    Vout = Jc_ini(:,feed_port) - X.' * (U.' * (M_q2ql.'* Jb));
else 
    Vout = Jc_ini(:,feed_port) - X.' * (U.' * Jb);
end


Vout = gather(Vout);