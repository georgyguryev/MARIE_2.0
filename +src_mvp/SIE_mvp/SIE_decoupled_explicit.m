function [Icu, Ics] = SIE_decoupled_explicit(Jb_cu, Jb_cs, b_Icu, X_cu, Z_L, rhs_cp)
% function [Vout] = SIE_decoupled_explicit(Jb, Jc_o, Zc, Zbc, feed_port)
% computes coil surface currents from inital guess (SIE solution) and
% contribution of tissue currents (dense coupling)
% ------------------------------------------------------------------------
% Input:
%        - Jb_cu - unscaled tissue polarization currents
%        - Jb_cs - scaled equivalent tissue currents
%        - b_Icu - port currents due to applied voltage in free space
%        - X_cu - operator that maps tissue currents to perturb. of port
%                  currents due to the presence of the body
%        - Z_L - matching load attached to the excitation ports
%        - rhs_cp - actual voltages applied to coil ports 
% Output:
%        - Icu - resulting port currents (w.o. matching loads)
%        - Ics - resulting port currents (w. matching loads)
% ------------------------------------------------------------------------ 
% Author: Georgy Guryev, Cambridge, MA, 2019    
 

Np = size(Jb_cu,2);
I  = eye(Np);

Icu = b_Icu + X_cu * Jb_cu;
Ics = (I - inv(I + inv(b_Icu * Z_L)))*(b_Icu * rhs_cp + X_cu * Jb_cs);

Icu = gather(Icu);
Ics = gather(Ics);