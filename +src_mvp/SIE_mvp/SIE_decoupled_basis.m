function [Icu, Ics]  = SIE_decoupled_basis(Jb_cu, Jb_cs, X_cu, b_Icu, Z_L, rhs_cp, M_q2ql)
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

Np = size(Jb_cu,2);
I  = eye(Np);

if nargin == 7  && size(X_cu,2) ~= size(M_q2ql,1)
    Icu = b_Icu + X_cu * (M_q2ql.' * Jb_cu);    
    Ics = (I - inv(I + inv(b_Icu * Z_L)))*(b_Icu * rhs_cp + X_cu * (M_q2ql.' * Jb_cs));
else 
    Icu = b_Icu + X_cu * Jb_cu;
    Ics = (I - inv(I + inv(b_Icu * Z_L)))*(b_Icu * rhs_cp + X_cu * Jb_cs);
end

Icu = gather(Icu);
Ics = gather(Ics);