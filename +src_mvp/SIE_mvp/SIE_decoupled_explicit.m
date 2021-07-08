function Vout = SIE_decoupled_explicit(Jb, Jc_o, Zc, Zbc, feed_port)
% function [Vout] = SIE_decoupled_explicit(Jb, Jc_o, Zc, Zbc, feed_port)
% computes coil surface currents from inital guess (SIE solution) and
% contribution of tissue currents (dense coupling)
% ------------------------------------------------------------------------
% Input:
%        - Jb - tissue polarization currents (vector)
%        - Jc_o - 
%        - Zc - coil-to-coil interaction matrix
%        - Zbc - coil-to-body interaction matrix
%        - feed_port - index of the feed/excitation port
% Output:
%        - Vout - the resulting solution for surface currents 
% ------------------------------------------------------------------------ 
% Author: Georgy Guryev, Cambridge, MA, 2019    

Vout = Jc_o(:,feed_port) - Zc \ (Zbc.' * Jb);

Vout = gather(Vout);