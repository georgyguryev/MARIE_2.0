function [Vout] = VIE_mvp(Jb, mvp, scatterer, dims, S_ql, freq)
% function [Vout] = VIE_mvp(Jb, mvp, scatterer, dims, S_ql, freq)
% implements the matrix-vector product for JVIE-I formulation
% (JVIE-2 when preconditioned) 
% ------------------------------------------------------------------------
% Input:
%        - Jb - tissue polarization currents (vector)
%        - mvp - instance of class with core operator mvps (Nop, Kop, etc.)
%        - scatterer - object that contains all tissue specific information
%        - dims - object keeps all dimensions of the problem
%        - S_ql - scatterer masks
%        - freq - working frequency 
% Output:
%        - Vout - result of the matrix-vector product with tissue currents 
% ------------------------------------------------------------------------ 
% Author: Georgy Guryev, Cambridge, MA, 2019    



% get electromagnetic properties
emu = src_utils.EM_utils(freq);

% get material properties
Mcr = scatterer.prop_vie.Mcr;
Mrc = scatterer.prop_vie.Mrc;

% preallocate tensor of body currents
Eout = zeros(dims.vie);
Jin  = zeros(dims.vie);

% copy scatterer polarization curents to tensor format
Jin(S_ql) = Jb;

%% apply N operator to input currents
E_b = mvp.N_vie(Jin);

% reshape Electric fields and currents for Hadamard product
E_b = reshape(E_b, dims.vie);
Jin = reshape(Jin, dims.vie);

% multiply J_tot by Gram matrix
G_J = mvp.G_vie(Jin);
G_J = reshape(G_J, dims.vie);

% multiply by material properties
for i = 1:dims.ql
    
    % compute the JVIE I formulation
    Eout(:,:,:,i) = (1 / emu.ce) .* Mrc .* (G_J(:,:,:,i) - Mcr .* E_b(:,:,:,i));
    
end

% select scatterer unknowns
Vout = Eout(S_ql);
