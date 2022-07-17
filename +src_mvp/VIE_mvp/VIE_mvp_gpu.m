function [Vout] = VIE_mvp_gpu(Jb, Jin, mvp, scatterer, dims, S_ql, freq)
% function [Vout] = VIE_mvp_gpu(Jb, mvp, scatterer, dims, S_ql, freq)
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

% transfer input vector to GPU
% Jb = gpuArray(Jb);

% get electromagnetic properties
emu = src_utils.EM_utils(freq);

% get material properties
Mcr = scatterer.prop_vie.Mcr;
Mrc = scatterer.prop_vie.Mrc;

% preallocate tensor of body currents
% Eout = gpuArray(zeros(dims.vie));
% Jin  = gpuArray(zeros(dims.vie));

% Eout = zeros(dims.vie);
% Jin = zeros(dims.vie);

% Eout = gpuArray(Eout);
% Jin = gpuArray(Jin);


% copy scatterer polarization curents to tensor format
Jin(S_ql) = Jb;

%% apply N operator to input currents

E_b = mvp.N_vie_gpu(Jin);

% reshape Electric fields and currents for Hadamard product
E_b = reshape(E_b, dims.vie);
Jin = reshape(Jin, dims.vie);

% multiply J_tot by Gram matrix
G_J = mvp.G_vie(Jin);
G_J = reshape(G_J, dims.vie);

% multiply by material properties
for i = 1:dims.ql
    
    % compute the JVIE II formulation
    E_b(:,:,:,i) =  1 / emu.ce .* Mrc .* (G_J(:,:,:,i) - Mcr .* E_b(:,:,:,i));

end

% select scatterer unknowns
% Vout = gather(E_b(S_ql));
Vout = E_b(S_ql);

%% deallocate memory 

% clear Eout Jin Jb E_b G_J Mcr Mrc;
