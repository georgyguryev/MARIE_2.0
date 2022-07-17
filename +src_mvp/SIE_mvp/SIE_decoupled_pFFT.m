function [Icu, Ics] = SIE_decoupled_pFFT(Jb_cu, Jb_cs, X_cu, b_Icu, L_ext, S, Z_L, rhs_cp)
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

% emu = src_utils.EM_utils(freq);



% L_b = L_ext(gpuArray(S * Jb));

Jb_cu_ext = S * Jb_cu;
Jb_cs_ext = S * Jb_cs;

N_feeds = size(Jb_cu_ext,2);

I = eye(N_feeds);
Ib_cu = zeros(size(Jb_cu_ext));
Ib_cs = zeros(size(Jb_cs_ext));

for port = 1:N_feeds
    Ib_cu(:,port) = gather(reshape(L_ext(Jb_cu_ext(:, port)), [],1));
    Ib_cs(:,port) = gather(reshape(L_ext(Jb_cs_ext(:, port)), [],1));
end

% Vout = Jc_o(:,feed_port) - 1 / emu.ce * Zc_inv * (P.' * L_b(:));
Icu = b_Icu + X_cu * Ib_cu;

Ics = (I - inv(I + inv(b_Icu * Z_L))) * (b_Icu * rhs_cp + X_cu * Ib_cs);

