function [E,H] = PlaneWave_Excitation(scatterer,k,omega,dims,Eo)

% PlaneWave_Excitation(r,k,omega_mu,Eo)
%%    Function to Generate a Plane wave excitation
% _________________________________________________________________________
%
%       Generates the E and H fields of a plane wave
%       The electric field magnitude is 1V/m
% _________________________________________________________________________
%
%% INPUT
%   r               4D (LxMxNx3) array with domain voxelized grid coordinates
%   k               vector with wavenumbers in the medium [kx ky kz]
%   omega_mu        omega*mu
%   Eo              polarization vector for E-field
%
%% OUTPUT
%   E               Electric field (LxMxNx3)
%   H               Magnetic field (LxMxNx3)
%
% -------------------------------------------------------------------------
%
%%   This function is part of MARIE
%   MARIE - Magnetic Resonance Integral Equation suite
%           Jorge Fernandez Villena   -- jvillena@mit.edu
%           Athanasios G. Polimeridis -- thanos_p@mit.edu
%           Copyright ï¿½ 2014
%           RLE Computational Prototyping Group, MIT
% 
%           This software is free and open source
%           Distributed under the GNU-GPLv3 terms
%           For details see MARIE_license.txt
%
% _________________________________________________________________________


% -------------------------------------------------------------------------
% Obtain voxel resolution and domain size
% -------------------------------------------------------------------------

res = scatterer.dom_vie.res;

x = scatterer.dom_vie.x_tensor;
y = scatterer.dom_vie.y_tensor;
z = scatterer.dom_vie.z_tensor;

L = dims.L_vie;
M = dims.M_vie;
N = dims.N_vie;

dx = x(2,1,1) - x(1,1,1);
dy = y(1,2,1) - y(1,1,1);
dz = z(1,1,2) - z(1,1,1);

% -------------------------------------------------------------------------
% polarization H-field
% -------------------------------------------------------------------------
mu_0 = 4 * pi * 1e-7;

Ho = 1 / (omega * mu_0) * cross(k,Eo);

% -------------------------------------------------------------------------
% extract wavenumbers
% -------------------------------------------------------------------------
kx = k(1);
ky = k(2);
kz = k(3);

% -------------------------------------------------------------------------
% traveling direction of excitation
% -------------------------------------------------------------------------
if kx ~= 0
    Ex =  exp(-1j * kx * x ) *  (exp(-1j*kx*dx/2) - exp(1j*kx*dx/2)) / (-1j*kx);
else
    Ex = dx * ones(L,M,N);
end

if ky ~= 0
    Ey = exp(-1j * ky * y ) *  (exp(-1j*ky*dy/2) - exp(1j*ky*dy/2)) / (-1j*ky);
else
    Ey = dy * ones(L,M,N);
end

if kz ~= 0
    Ez = exp(-1j * kz * z ) *  (exp(-1j*kz*dz/2) - exp(1j*kz*dz/2)) / (-1j*kz);
else
    Ez = dz * ones(L,M,N);
end
%
Eexc = Ex .* Ey .* Ez;

% -------------------------------------------------------------------------
% apply polarization
% -------------------------------------------------------------------------

E = [Eo(1)*Eexc(:) ; Eo(2)*Eexc(:) ; Eo(3)*Eexc(:)];
%
H = [Ho(1)*Eexc(:) ; Ho(2)*Eexc(:) ; Ho(3)*Eexc(:)];

% -------------------------------------------------------------------------
% reshape and scale to 1V/m
% -------------------------------------------------------------------------

E = reshape(E, L, M, N, 3);
H = reshape(H, L, M, N, 3);


