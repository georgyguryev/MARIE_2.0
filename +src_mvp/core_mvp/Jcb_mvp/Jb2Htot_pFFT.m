function H_tot = Jb2Htot_pFFT(Jb, K_ext, L_ext, Zc2b_Kop, Zc2b_T, Zc_inv_hat, P,S, b_K, freq)
% function Jb2Htot() maps equivalent tissue currents Jb to total magnetic fields within the 
% tissue domain
%
% Input:
%       Jb - equivalent tissue currents
%       K_ext - magnetic K operator defined over the extended domain 
%       L_ext - electric L operator over the extended domain (coil + body)
%       Zc2b_Kop - precorrection magnetic coupling matrix 
%       Zc2b_T   - precorrection electric coupling matrix
%       Zc_inv_hat - inverted coil-to-coil interaction matrix (with ports
%                    loaded with the external matching circuits)
%       P - projection matrix that maps coil surface currents to equivalent
%           volumetric currents defined in the extended domain
%       S - projection matrix that maps volumetric tissue currents to the
%           volumetric currents in the extended domain (simple indexing)
%       b_K - precomputed expression defined as b_K = Zc_inv_hat * F * rhs_cp;
%       freq - working frequency 

% get EM parameters
emu = src_utils.EM_utils(freq);


Nfeed = size(b_K,2);
H_tot = zeros(size(Jb));

for port = 1:Nfeed
    L_b = L_ext(gpuArray(S * Jb(:,port)));
    
    % compute scaled surface currents 
    Jc  = b_K(:, port) -  Zc_inv_hat * (1 / emu.ce * P.' * L_b(:) + Zc2b_T * Jb(:,port));
    
    % find total equivalent currents 
    Jtot      = P * Jc + S * Jb(:,port);
    
    % compute total magnetic field within the tissue before precorrection
    H_vox_ext = K_ext(Jtot);
    
    % precorrect the total magnetic field
    H_tot(:,port) = gather(Zc2b_Kop * Jc +  S.' * H_vox_ext(:));
end
 