function Vout = VSIE_coupled_pFFT(Jcb, mvp, scatterer, dims, precorrection, P)

%% separate coil from body currents

Jc = Jcb(1:dims.N_sie, 1);
Jb = Jcb(dims.N_sie + 1:end, 1);
E  = zeros(size(Jcb,1),3);

% construct Jcb with zero coil currents
Jco_b = [zeros(size(Jc)); Jb];
 
% map Jb to Jb_ext
Jb_ext = mvp.Jcb2Jtot(Jco_b);

%% Alternative with Zc_vox precorrection

% project currents on extended voxelized domain
J_tot = mvp.Jcb2Jtot(Jcb);

% apply L operator to total currents
L_b = mvp.L_ext(J_tot);

% reshape tensor to the vector dimensions
L_b = reshape(L_b, size(J_tot));

% [Ecb, Ebc] = fJtot2Jcb(L_b);
E(:,1) = mvp.Jtot2Jcb(L_b);

%% fill Vout_b

Jb_ext = reshape(Jb_ext, dims.ext);

% compute scaling factor for tissue currents
Mc_inv_Gram = scatterer.prop_ext.Mc_inv * scatterer.dom_ext.res.^3;

% scale Jb_ext by material properties
Ebb_4d = Mc_inv_Gram .* Jb_ext;

% interpolate Jb_ext on body and coil domains
E(:,2) = mvp.Jtot2Jcb(Ebb_4d(:));

% add precorrection term
E(:,3) = precorrection.Z_tot * Jcb;

% reshape resulting expression and scale by P
E = reshape(E, 3 * size(Jcb,1),1);
Vout = P * E;
