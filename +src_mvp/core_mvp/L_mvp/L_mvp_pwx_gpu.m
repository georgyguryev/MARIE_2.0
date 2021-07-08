function Jout = L_mvp_pwx_gpu(Jin, fN_mvp, fG_mvp, dom_dims)
% function L_mvp_pwx(Jin, fN, fG, dom_dims)

%% apply Nop:  N(Jin)

Jout = fN_mvp(Jin);

% reshape Electric fields and currents for Hadamard product
Jout = reshape(Jout, dom_dims);
Jin  = reshape(Jin, dom_dims);

% multiply J_tot by Gram matrix
G_J = fG_mvp(Jin);
G_J = reshape(G_J, dom_dims);

Jout = Jout - G_J;

% gather reluts
clear Jin;
