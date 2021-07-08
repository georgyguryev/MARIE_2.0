function [Hinc] = Jc2K_coupling(Jc,r,SCOIL, pwx, freq)
%%    Compute fields due to the solution to the SIE problem
% _________________________________________________________________________
%
%   Applies operators to generate the fields and other figures of merit
%
% _________________________________________________________________________
%
%% INPUT
%       SCOIL structure
%           index - mapping of the internal edge number to dof number
%           etod - etod numbering of elements and edges
%           node - coordinates of the nodes 
%           edge - numbering of edges
%           elem - 3 indexes of the nodes defining an element
%           Index_elem - mapping of index to elements
%           port - port definition
%           ...
%       Jc - currents in the coil
%       freq - frequency
%       r - domain
%
%
%% OUTPUT
%       H        Solution magnetic field (LxMxNx3)
%
%
% -------------------------------------------------------------------------
%
%   J. Fernandez Villena -- jvillena@mit.edu
%   A.G. Polimeridis -- thanos_p@mit.edu
%   Computational Prototyping Group, RLE at MIT
%
% _________________________________________________________________________




% -------------------------------------------------------------------------
%                 define EM vars and constants
% -------------------------------------------------------------------------

[L,M,N,~] = size(r);

% x,y,z coordinates of the domain
xd = r(:,:,:,1);
yd = r(:,:,:,2);
zd = r(:,:,:,3);

% Coordinates of all the domain
Dcoord = [xd(:), yd(:), zd(:)];

% -------------------------------------------------------------------------
%         Compute incident fields due to Jc
% -------------------------------------------------------------------------


[Hinc] = src_coupling.SVIE_M_Coil2Scat_PM_par(Dcoord,SCOIL.index,SCOIL.etod,SCOIL.Ct,SCOIL.Ln,SCOIL.Pn,freq,Jc);

res = abs(xd(2,1,1) - xd(1,1,1));
Hinc = res^3 * Hinc;

Hinc = reshape(Hinc,L,M,N,3);
    
if pwx == 12
   Hinc_temp = zeros(L,M,N,12);
   Hinc_temp(:,:,:,1) = Hinc(:,:,:,1);
   Hinc_temp(:,:,:,5) = Hinc(:,:,:,2);
   Hinc_temp(:,:,:,9) = Hinc(:,:,:,3);
   Hinc = Hinc_temp;
end
% -------------------------------------------------------------------------
%                 And it is done
% -------------------------------------------------------------------------


