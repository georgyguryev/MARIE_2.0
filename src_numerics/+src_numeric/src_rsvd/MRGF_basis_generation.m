function  [Uin,Xin,Pin,Dcoord] = MRGF_basis_generation(RHBM,freq,idxI,Lmax,tol,Uin_full_fname, Uin_fname, PXin_fname,GPU_flag,deim_type)
%
% _________________________________________________________________________
%
%   Function that generates the basis and DEIM points of MRGFs
%   for a given RHBM and idxI, for a number of excitations
%
%   INPUT:  RHBM            realistic human body model struct
%           freq            frequency
%           idxI            coil domain: indexes of the dipole positions in domain
%           Lmax            number of vectors used as excitations
%           tol             relative tolerance for the basis truncation
%           Uin_full_fname  name of temp file with Uin bassis before trunc.
%           Uin_fname       name of temp file with Uin after truncation
%           PXin_fname      name of temp file with DEIM Xin, Pin matrices
%                           (see deim_type) 
%           GPU_flag        1 to use the GPU
%           deim_type       0 - if DEIM is aplied  before truncation of Uin, 
%                           1 - if DEIM is aplied  after truncation of Uin 
%
%   OUTPUT: Uin       basis vectors of the indicent fields
%           Xin       deim matrix
%           Pin       incidence matrix of deim points 
%           Dcoord    coordinates of deim points
%           
% _________________________________________________________________________
%
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
   Uin_fname = [];
   PXin_fname = [];
   GPU_flag = 0;
   deim_type = 0;
end
if(nargin < 6 )
%    temp_prefix = [ ];
   Uin_fname = [];
   PXin_fname = [];
   deim_type = 0;
   GPU_flag = 0;
end
if(nargin < 7)
    PXin_fname = [];
    GPU_flag = 0;
    deim_type = 0;
end;
if(nargin < 8 || isempty(GPU_flag))
    GPU_flag = 0;
    deim_type = 0;
end

% -------------------------------------------------------------------------
%                   Define local file names 
% -------------------------------------------------------------------------

% if ~isempty(temp_prefix)
%     Uin_fname         = sprintf('%s_MRGF_Uin_tol_%1.e.mat',temp_prefix, tol);
%     PXin_before_fname = sprintf('%s_MRGF_DEIM_before_tr_tol_%1.e.mat',temp_prefix, tol);
%     PXin_after_fname  = sprintf('%s_MRGF_DEIM_after_tr_tol_%1.e.mat',temp_prefix, tol);
% end;
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


% -------------------------------------------------------------------------
% Initialization of variables
% -------------------------------------------------------------------------

fid = 1;

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
% Run Randomized SVD for basis construction (if it does not exist)
% -------------------------------------------------------------------------

if (exist(Uin_full_fname,'file') && ~isempty(Uin_full_fname)...
        && exist(Uin_full_fname,'file') && ~isempty(Uin_full_fname))
    
    if(0 == deim_type)
        % load existing full basis Uin_full
        load(Uin_full_fname);
    else
        %load existing truncated basis Uin
        load(Uin_fname);
    end;
else
    % compute approximate numerical basis for left subspace of coupling matrix
    [Uin,~,~,numsv] = rSVD_MRGF(Fdirect,Fadjoint,RHBM.r,idxI3,idxS3,Lmax,tol,blocksize);
    save(Uin_full_fname, 'Uin', 'numsv','-v7.3');
end;
%%
% -------------------------------------------------------------------------
% Apply the DEIM procedure on the incident basis
% -------------------------------------------------------------------------

if (0 == deim_type)
    tinter = tic;
    [Pin,Xin,Dcoord] = deim_MRGF(Uin,numsv,r,idxS);
    timeinc = toc(tinter);
    
    %print DEIM info
    ndeim = size(Pin,2);
    fprintf(fid, '\n\n DEIM algorithm applied on %d vectors. Elapsed time %g \n', ndeim, timeinc);
    
end;

% -------------------------------------------------------------------------
% Truncate resulting Basis and interpolation matricies
% -------------------------------------------------------------------------

% do the actual truncation of the basis
Uin = Uin(:,1:numsv);

if(0 == deim_type)
    % truncate also the DEIM Xdx extended matrix
    Xin = Xin(1:numsv,:);
end;

% -------------------------------------------------------------------------
% Save truncated Basis Uin and matricies Xin, Pin
% -------------------------------------------------------------------------

if ~isempty(Uin_fname)
    % save truncated basis
    save(Uin_fname, 'Uin', 'numsv','-v7.3');
end;

% -------------------------------------------------------------------------
% Apply DEIM to truncated basis Uin; 
% -------------------------------------------------------------------------

% if user selected this option 
% apply DEIM, save results to temp file 
if (1 == deim_type)
    
    tinter = tic;
    [Pin,Xin,Dcoord] = deim_MRGF(Uin,numsv,r,idxS);
    timeinc = toc(tinter);
    
    % print DEIM info
    ndeim = size(Pin,2);
    fprintf(fid, '\n\n DEIM algorithm applied on %d vectors. Elapsed time %g \n', ndeim, timeinc);  
end;

% Save matricies Xin, Pin
if ~isempty(PXin_fname)
    save(PXin_fname, 'Pin','Xin', 'Dcoord','-v7.3');
end;
