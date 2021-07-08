function Escat = Jb2Escat_mvp(Jb, mvp, dims, freq)


%% compute total current and map volumetric current

% get elecromagnetic constants
emu = src_utils.EM_utils(freq);


% form a vector of body currents
Jcb = [zeros(dims.N_sie, 1); Jb];

% map Jc, Jb to J_tot
J_tot = mvp.Jcb2Jtot(Jcb);

%% compute L operator for L(Jb_ext) and for L(Jtot)

% Compute fields, produced by coil surface currents within extended domain
Vb_vox = mvp.L_ext(J_tot);

% get fields within a scatterer domain
Vcb = mvp.Jtot2Jcb(Vb_vox(:));

Escat = 1 / emu.ce * Vcb(dims.N_sie+1:end);

gather (Escat);
    
