function Escat = Jb2Escat_vie_mvp(Jb, mvp, freq)

%% compute total current and map volumetric current

% get elecromagnetic constants
emu = src_utils.EM_utils(freq);

%% compute L operator for L(Jb_ext) and for L(Jtot)

% Compute fields, produced by coil surface currents within extended domain
Vb2b = mvp.L_vie(Jb);

Escat = 1 / emu.ce * Vb2b;

gather (Escat);
    
