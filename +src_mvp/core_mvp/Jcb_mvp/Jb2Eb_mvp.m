function Eb = Jb2Eb_mvp(Jb, dims, scatterer, freq)

%% compute total current and map volumetric current

% get elecromagnetic constants
emu = src_utils.EM_utils(freq);

%% compute L operator for L(Jb_ext) and for L(Jtot)

% get fields within a scatterer domain
Mc_inv = scatterer.prop_vie.Mc_inv;
idxS = scatterer.index_vie.index_ql(dims.ql, dims.Nvox_vie);

Jb_tensor = zeros(dims.vie);
Jb_tensor(idxS) = Jb;

Eb = zeros(size(Jb_tensor));

for i =1:dims.ql

    Eb(:,:,:,i) = 1 / emu.ce * Mc_inv .* Jb_tensor(:,:,:,i);    
end
