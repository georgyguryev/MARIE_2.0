function [Vout] =  VSIE_decoupled_pFFT(Jb, L_ext, dims, Mc_inv_Gram, Zbc_prec, Zbc_prec_T, Zc_inv, P,S, freq)

%% separate coil from body currents
emu = src_utils.EM_utils(freq);

Jb = gpuArray(Jb); 

% map Jb to Jb_ext
Jb_ext = S * Jb;
Jb_ext = gpuArray(Jb_ext);

%% Alternative with Zc_vox precorrection

% apply L operator to total currents
L_b = L_ext(Jb_ext);

% reshape tensor to the vector dimensions
L_b = L_b(:);

Jb_ext = reshape(Jb_ext, dims.ext);

% scale Jb_ext by material properties
E_0 = 1 / emu.ce * Mc_inv_Gram .* Jb_ext; 

S_0 = S.' * E_0(:);
S_1 = - 1 / emu.ce * S.' * L_b;
X_0 = Zbc_prec_T * Jb + 1 / emu.ce * P.' * L_b;
X_1 = Zc_inv * X_0;
X_2 = L_ext(P * X_1);
S_23 = Zbc_prec * X_1 + 1 / emu.ce * S.' * X_2(:);

Vout = S_0 + S_1 + S_23;
Vout = gather(Vout);