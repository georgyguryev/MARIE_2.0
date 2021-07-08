function  [Uin,Xin,Pin,Dcoord] = MRGF_basis_generation_after_tr(RHBM,freq,idxI,Lmax,tol,temp_name,GPU_flag)
%
% _________________________________________________________________________
% _________________________________________________________________________
%
%   Function that generates the basis and DEIM points of MRGFs
%   for a given RHBM and idxI, for a number of excitations
%
%   INPUT:  RHBM      realistic human body model struct
%           freq      frequency
%           idxI      coil domain: indexes of the dipole positions in domain
%           Nexc      number of vectors used as excitations
%           tol       relative tolerance for the basis truncation
%           temp_name name for the file of temporary storage of elements
%           GPU_flag  1 to use the GPU
%
%   OUTPUT: Uin       basis vectors of the indicent fields
%           Xin       deim matrix
%           Pin       incidence matrix of deim points 
%           Dcoord    coordinates of deim points
%
% _________________________________________________________________________
%
%
% -------------------------------------------------------------------------
%
%%   This function is part of MARIE
%   MARIE - Magnetic Resonance Integral Equation suite
%           Jorge Fernandez Villena   -- jvillena@mit.edu
%           Athanasios G. Polimeridis -- thanos_p@mit.edu
%           Copyright © 2016
%           RLE Computational Prototyping Group, MIT
% 
%           This software is free and open source
%           Distributed under the GNU-GPLv3 terms
%           For details see MARIE_license.txt
%
% _________________________________________________________________________

% Initialize variables.
if(nargin < 4 )
   fprintf(1, '\n ERROR: not enough arguments\n');
   return
end
if(nargin < 5 || isempty(tol))
   tol = 1e-4;
end
if(nargin < 6 )
   temp_name = [ ];
end
if(nargin < 7 || isempty(GPU_flag))
   GPU_flag = 1;
end


% -------------------------------------------------------------------------
%                   initialize stuff
% -------------------------------------------------------------------------

% generate EM constants
EMconstants;

% coordinates of domain
xd = RHBM.r(:,:,:,1);
yd = RHBM.r(:,:,:,2);
zd = RHBM.r(:,:,:,3);
[L,M,N,~] = size(RHBM.r);
nD = L*M*N; % number of voxels in the complete domain

% scatterer
idxS = find(abs(RHBM.epsilon_r(:)-1 + RHBM.sigma_e(:)) > 1e-12); % non-air positions
idxS3 = [idxS; nD+idxS; 2*nD+idxS]; % the vector of non-air positions in 3D grid
nS = length(idxS);
xs = xd(idxS);
ys = yd(idxS);
zs = zd(idxS);
Scoord = [xs(:), ys(:), zs(:)];


% dipoles
if isempty(idxI)
    idxI = find((abs(e_r(:))-1 )== 0); % get indexes of air elements
end
idxI3 = [idxI; nD+idxI; 2*nD+idxI]; % the vector of input positions in 3D grid
nI = length(idxI);
xi = xd(idxI);
yi = yd(idxI);
zi = zd(idxI);
Icoord = [xi(:), yi(:), zi(:)];


% -------------------------------------------------------------------------
%                 define EM vars and constants
% -------------------------------------------------------------------------

% properties at the given frequency
e_r = RHBM.epsilon_r - 1j*RHBM.sigma_e/(eo*omega);

% to simplify
r = RHBM.r;

% get the voxel side
dx = r(2,1,1,1) - r(1,1,1,1);
% Gram matrix (value)
Gram = dx^3; % volume of the voxel

% -------------------------------------------------------------------------
%                  Generate circulants
% -------------------------------------------------------------------------

% compute the circulants
[RHBM.fN] = getOPERATORS(r,freq,'N',[],'DEMCEM');
%[RHBM.fK] = getOPERATORS(r,freq,'K',[],'DEMCEM');


% -------------------------------------------------------------------------
% Initialization of variables
% -------------------------------------------------------------------------

fid = 1;
% temporary storage in  .\body_data folder
if ispc
    outpath = strcat('.\src_mrgf\test_generated_basis\', temp_name);
else
    outpath = strcat('./src_mrgf/test_generated_basis/', temp_name);
end

% -------------------------------------------------------------------------
% Prepare handles and supplementary parameters for rSVD
% -------------------------------------------------------------------------

% define increment block size
blocksize = 50;

Jvec = zeros(L,M,N,3);

% create handles to compute produced fields Einc
Fdirect  = @(J) E_field_Nop(J, RHBM.fN,Gram,freq,[],GPU_flag);
Fadjoint = @(J) E_field_Nop(conj(J),RHBM.fN,Gram,freq,[],GPU_flag);
% -------------------------------------------------------------------------
% Run Randomized SVD for basis construction
% -------------------------------------------------------------------------

% generate Q - approx of operator's left subspace

% [Q,~] = getQ_MRGF(Fdirect,RHBM.r,3*nS,3*nI,idxI3,idxS3,Lmax,tol,blocksize);

% [W,~] = getQ_MRGF(Fadjoint,RHBM.r,3*nI,3*nS,idxS3,idxI3,Lmax,tol,blocksize);


% compute approximate numerical basis for left subspace of coupling matrix
[Uin,~,~,numsv] = rSVD_MRGF(Fdirect,Fadjoint,RHBM.r,idxI3,idxS3,Lmax,tol,blocksize);


if (~isempty(temp_name))
    Basis_file = sprintf('%s_MRGF_Uin.mat', outpath);
    save(Basis_file, 'Uin', '-v7.3');
end


%%



% -------------------------------------------------------------------------
% Apply the DEIM procedure on the incident basis
% -------------------------------------------------------------------------

tinter = tic;

[Pin,Xin,Dcoord] = deim_MRGF(Uin,numsv,r,idxS);

timeinc = toc(tinter);

ndeim = size(Pin,2);
fprintf(fid, '\n\n DEIM algorithm applied on %d vectors. Elapsed time %g \n', ndeim, timeinc);


% -------------------------------------------------------------------------
% Truncate resulting Basis and interpolation matricies
% -------------------------------------------------------------------------

% do the actual truncation of the basis
Uin = Uin(:,1:numsv);

% truncate also the DEIM Xdx extended matrix
Xin = Xin(1:numsv,:);


if (~isempty(temp_name))
    
    DEIM_file = sprintf('%s_MRGF_DEIM.mat', outpath);
    save(DEIM_file, 'Xin', 'Pin', 'Dcoord', '-v7.3');

    IBasis_file = sprintf('%s_MRGF_Utrunc.mat', outpath);
    save(IBasis_file, 'Uin', '-v7.3');
    
end

%% recompute Zin and Pin for the truncated basis Uin

clear Pin Xin Dcoord;

tinter = tic;

[Pin,Xin,Dcoord] = deim_MRGF(Uin,numsv,r,idxS);

timeinc = toc(tinter);

if (~isempty(temp_name))
    
    DEIM_file = sprintf('%s_MRGF_DEIM_after_tr.mat', outpath);
    save(DEIM_file, 'Xin', 'Pin', 'Dcoord', '-v7.3');    
end
