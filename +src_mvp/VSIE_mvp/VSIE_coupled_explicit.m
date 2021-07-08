function Vout = VSIE_coupled_explicit(Jcb, mvp, Zc, Zbc, Zbc_T)

%% separate coil from body currents

Jc = Jcb(1:mvp.dims.N_sie, 1);
Jb = Jcb(mvp.dims.N_sie + 1:end, 1);

%% compute C2C interactions

E_c2c = Zc * Jc;

%% compute UWU.' interactions

E_b2c = Zbc_T * Jb; 
E_c2b = Zbc * Jc;

%% compute B2B interactions

Jb    = gpuArray(Jb);
E_b2b = mvp.VIE(Jb);

%% compute the resulting field variation 
clear Jb;

Vout = [E_c2c  + E_b2c; -E_c2b + E_b2b];
Vout = gather(Vout);
