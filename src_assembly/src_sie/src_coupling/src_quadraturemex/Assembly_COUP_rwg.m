function [Zbc] = Assembly_COUP_rwg(Scoord,coil,dims,task_settings,freq,res,i_sie, idx_imp_near)
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
%       Zbc - Matrix (No x 3) with the contribution of each edge of the element
%
%
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
%            Define EM constants
% -------------------------------------------------------------------------


% get coil-related data
etod  = coil.etod;
node  = coil.node;
elem  = coil.elem;
index = coil.index;

% % list of near interactions
GL_order_sie  = task_settings.vsie.Np_quad_coup_sie;
GL_order_vie  = task_settings.vsie.Np_quad_coup_vie;

%% setup electomagnetic parameters
emu = src_utils.EM_utils(freq);
% -------------------------------------------------------------------------
% Define variables and allocate space
% -------------------------------------------------------------------------

RO = Scoord.';


% find physical edge for given dof(i_sie)
phys_edge = find(index == i_sie);

% find tirangles that share this physical edge
[idx_tr1,tr_1] = find(etod == phys_edge);
[idx_tr2,tr_2] = find(etod == -phys_edge);

% set of triangles that share common edge
tr = [tr_1, tr_2];

% extract verticies of supporting triangles
R1 = node(:,elem(1,tr)); % 3xNe with coordinates of the first node of all elements
R2 = node(:,elem(2,tr)); % 3xNe with coordinates of the second node of all elements
R3 = node(:,elem(3,tr)); % 3xNe with coordinates of the third node of all elements

%% compute length of common edge

% find nodes that form common edge
edge_node = setdiff([1,2,3], idx_tr1);

L1 = norm(node(:,elem(edge_node(2),tr_1)) - ...
          node(:,elem(edge_node(1),tr_1)), 2);

%%
    
% get the sign and index for the contribution
% MULT = [1,-1]; % +1 or -1

% define dimensions for coupling matrix assembly
NO = size(Scoord,1); % number of observation points
NE = size(tr,2); % number of elements
NC = size(i_sie,1); % number of dofs

% setup edge multiplier
MULT = etod(:,tr) ./ abs(etod(:,tr));

% setup column shift in coupling matrix
IDX = zeros(size(R1));
IDX(idx_tr1, 1) = 1;
IDX(idx_tr2, 2) = 1;

% -------------------------------------------------------------------------
% 1D cubature's number of points
% -------------------------------------------------------------------------

% LEVEL_DVrule = 5; % usually more than enough
[ Np_2D, Z1, Z2, Z3, wp ] = dunavant_rule (GL_order_sie);

% [Np_2D,wp,Z1,Z2,Z3] = Gauss_2Dt(1);
% -------------------------------------------------------------------------
% loop on the elements and fill the matrix
% -------------------------------------------------------------------------
%% Uncomment omp implementation when use on cluster!!!

if (dims.l == 4)
    % matlab implementation of linear coupling 
    [Zbc] = Assemble_coup_lin(R1,R2,R3,RO,NO,IDX,L1,emu.k0,res,Np_2D,Z1,Z2,Z3,wp,GL_VIE_order);
    scalefactor =  1 / emu.ce;
else
% call the mex fun: already openmpi-ed if available
% if (exist('ompQuadCoil2Scat', 'file') == 3)
%     [Zbc] = ompQuadCoil2Scat(R1(:),R2(:),R3(:),NE,RO(:),NO,IDX(:),MULT(:),NC,ko,Np_2D,Z1,Z2,Z3,wp);
% else

    if strcmp('struct', class(GL_order_vie))
        % compute matlab-based coupling quadrature assembly for constant
        % basis; here we use different order of quadratures to capture
        % field variation
        [Zbc] = Assemble_coup_const_imp(R1,R2,R3,RO,NO,IDX,L1,ko,res,Np_2D,Z1,Z2,Z3,wp,GL_order_vie, idx_imp_near);

          
%         keyboard
    else
        % compute matlab-based coupling quadrature assembly for constant basis
%         [Zbc] = Assemble_coup_const(R1,R2,R3,RO,NO,IDX,L1,ko,res,Np_2D,Z1,Z2,Z3,wp,GL_VIE_order);
        
        
%         % parallel compute matlab-based coupling quadrature assembly for constant basis
%         [Zbc] = Assemble_coup_const_par(R1,R2,R3,RO,NO,IDX,L1,ko,res,Np_2D,Z1,Z2,Z3,wp,GL_VIE_order);
    end
    
    % set a scaling factor for 
%     scalefactor =  1 / ce;


%     % mex-based quadrature function 
    [Zbc] = mexQuadCoil2Scat(R1(:),R2(:),R3(:),NE,RO(:),NO,IDX(:),MULT(:),NC,emu.k0,Np_2D,Z1,Z2,Z3,wp);
    
    scalefactor = - 1 / emu.ce / (4*pi);
end

%% scale resulting coupling matrix
% -------------------------------------------------------------------------
%             Final Z (with mult. constant) 
%
%      4*pi comes from omitting it in Green function!!
%       Z = discretization of  -e^scattered
% -------------------------------------------------------------------------

Zbc = scalefactor * Zbc;
