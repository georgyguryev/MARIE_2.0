function [Zbc] = Assembly_SCOUP_QMEX_par(Scoord,index,etod,node,elem,freq,LEVEL_DVrule) 
%%    Quadrature coupling for the SIE+VIE solver (PARALLEL Mexed)
% _________________________________________________________________________
%
%   Fucntion to generate the Quadrature Coupling SIE to VIE
%   Applies a parallel (MultiCore) mex function to compute all elements
%
% _________________________________________________________________________
%
%% Input
%       Scoord - coordinates of the observation points (No x 3)
%       node - coordinates of the nodes 
%       elem - 3 indexes of the nodes defining an element
%       etod - etod
%       index - mapping of the internal edge number to dof number
%       freq - frequency
%       LEVEL_DVrule - level for the DUNAVANT rule
%
%
%% Output
%       Zbc - Tensor (No x 3 x Nd) with the contribution of each edge of the element
%              Zbc(1000,3,200) is z contribution of 200-th edge to 1000-th element
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
%            Define EM constants
% -------------------------------------------------------------------------

mu = 4*pi*1e-7;
co = 299792458;
eo = 1/co^2/mu;
%
omega = 2 * pi * freq;
lambda  = co/freq;
ko = 2*pi/lambda;

% Free-space impedance
% eta = omega*mu/ko; %3.767303134617706e+002; 

% -------------------------------------------------------------------------
% Define variables and allocate space
% -------------------------------------------------------------------------

NO = size(Scoord,1); % number of observation points
NE = size(elem,2); % number of elements
NC = max(index); % number of dofs

RO = Scoord.';

R1 = node(:,elem(1,:)); % 3xNe with coordinates of the first node of all elements
R2 = node(:,elem(2,:)); % 3xNe with coordinates of the first node of all elements
R3 = node(:,elem(3,:)); % 3xNe with coordinates of the first node of all elements
    
% get the sign and index for the contribution
ABSNUM = abs(etod(:,:)); % internal index of the edge
MULT = etod(:,:)./ABSNUM; % +1 or -1
IDX = index(ABSNUM);


% -------------------------------------------------------------------------
% 1D cubature's number of points
% -------------------------------------------------------------------------

% LEVEL_DVrule = 5; % usually more than enough
[Np_2D, Z1, Z2, Z3, wp] = dunavant_rule ( LEVEL_DVrule );

% [Np_2D,wp,Z1,Z2,Z3] = Gauss_2Dt(1);
% -------------------------------------------------------------------------
% loop on the elements and fill the matrix
% -------------------------------------------------------------------------

% call the mex fun: already openmpi-ed if available
if (exist('ompQuadCoil2Scat', 'file') == 3)
    [Zbc] = ompQuadCoil2Scat(R1(:),R2(:),R3(:),NE,RO(:),NO,IDX(:),MULT(:),NC,ko,Np_2D,Z1,Z2,Z3,wp);
else
    [Zbc] = mexQuadCoil2Scat(R1(:),R2(:),R3(:),NE,RO(:),NO,IDX(:),MULT(:),NC,ko,Np_2D,Z1,Z2,Z3,wp);
end
Zbc = reshape(Zbc,NO,3,NC);

% Zbc = Zbc/(1j*ko);
% scalefactor = (1j*eta)/(4*pi*ko);
% -------------------------------------------------------------------------
%             Final Z (with mult. constant) 
%
%      4*pi comes from omitting it in Green function!!
%       Z = discretization of  -e^scattered
% -------------------------------------------------------------------------
ce = 1i*omega*eo;
scalefactor = - 1 / ce / (4*pi);

Zbc = scalefactor * Zbc;
    