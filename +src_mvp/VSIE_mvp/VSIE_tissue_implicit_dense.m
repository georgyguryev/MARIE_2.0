function [Vout] = VSIE_tissue_implicit_dense(Jc, mvp, vie_solver)
% Vout = VSIE_decoupled_implicit(Jc);
% function applies VSIE operator to coil current

% apply coupling operator to coil currents
Vbc = mvp.c2b_implicit_coupling(Jc);

% solve VIE problem
% Jb = vie_solver.run(Vbc,[]);


[Jb,rel_res,res_vec] = vie_solver.run(Vbc,[]);

Jb = gather(Jb);

% apply coupling operator to body currents+
Vcb = mvp.b2c_implicit_coupling(Jb);

%% compute MVP for nested system (coil excitation for coupled problem)

Vout = mvp.Z_coil * Jc + Vcb;