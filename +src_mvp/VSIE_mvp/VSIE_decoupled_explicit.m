function Vout = VSIE_decoupled_explicit(Jb, mvp, Zc_inv, Zbc)


% Jb = gpuArray(Jb);

%% compute field contribution of pertrubations in surface currents  

% E_b2b_pert = Zbc * (Zc \ (Zbc_T * Jb));

E_b2b_pert = Zbc * (Zc_inv * (Zbc.' * Jb));

%% compute B2B interactions

Jb    = gpuArray(Jb);
E_b2b = mvp.VIE(Jb);

%% compute the resulting field variation 

Vout =  E_b2b + E_b2b_pert;
Vout = gather(Vout);