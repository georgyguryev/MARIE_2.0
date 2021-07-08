function Vout = VSIE_decoupled_basis(Jb, mvp)

% separate coil from body currents
Jb = gpuArray(Jb);

% compute UWU.' interactions
E_UWU = mvp.uwu_implicit_coupling(Jb);

% compute B2B interactions
E_b2b = mvp.VIE(Jb);

% compute the resulting field variation 
Vout = gather(E_b2b + E_UWU);
% Vout = gather(E_b2b);

clear Jb E_UWU E_b2b;
