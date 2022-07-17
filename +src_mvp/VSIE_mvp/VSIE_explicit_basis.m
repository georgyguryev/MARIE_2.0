function Vout = VSIE_explicit_basis(Jcb, mvp, dims)


%% separate coil from body currents

% Jc = gpuArray(Jcb(1:dims.N_sie, 1));
% Jb = gpuArray(Jcb(dims.N_sie + 1:end, 1));

Jc = Jcb(1:dims.N_sie, 1);
Jb = Jcb(dims.N_sie + 1:end, 1);

%% compute C2C interactions

E_c2c = mvp.Z_coil * Jc;

%% compute C2B interactions

% E_c2b = 1 / emu.ce * Gram * U * (alpha * Jc);
E_c2b = mvp.c2b_implicit_coupling(Jc);

%% compute B2C interactions

% E_b2c = 1 / emu.ce * Gram * alpha.' * (U.' * Jb);
E_b2c = mvp.b2c_implicit_coupling(Jb);

%% compute B2B interactions

% E_b2b = mvp.VIE(Jb);

E_b2b = mvp.VIE(Jb);


%% compute the resulting field variation 

Vout = gather([E_c2c  + E_b2c; -E_c2b + E_b2b]);

clear Jc Jb E_b2c E_c2b E_b2b;
