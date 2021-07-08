function Vout = SIE_decoupled_pFFT(Jb, feed_port, Jc_o, Zc_inv, L_ext, P, S, freq)
% function [Vout] = SIE_decoupled_basis(Jb, feed_port, Jc_o, Zc_inv, L_ext, P, S, freq)
% computes coil surface currents from inital guess (SIE solution) and
% contribution of tissue currents (pFFT coupling)
% ------------------------------------------------------------------------
% Input:
%        - Jb - tissue polarization currents (vector)
%        - feed_port - index of the feed/excitation port
%        - Jc_o - initial guess for surface currents 
%        - Zc_inv - inverse  coil-to-coil interaction matrix
%        - L_ext - L opertor mvp defined whitin extended domain
%        - P,S  - coil and tissue projection matrices (to ext. domain)
%        - freq - working frequency
% Output:
%        - Vout - the resulting solution for surface currents 
% ------------------------------------------------------------------------ 
% Author: Georgy Guryev, Cambridge, MA, 2019 

emu = src_utils.EM_utils(freq);

L_b = L_ext(gpuArray(S * Jb));

Vout = Jc_o(:,feed_port) - 1 / emu.ce * Zc_inv * (P.' * L_b(:));

Vout = gather(Vout);