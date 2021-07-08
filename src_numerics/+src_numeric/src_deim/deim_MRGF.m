function [Pin,Xin,Dcoord] = deim_MRGF(Uin,numsv,r,idxS)
%% deim_MRGF applies Discrete empirical interpolation method

%   INPUT:
%
%
%   OUTPUT:
%   Pin  - DEIM selection matrix
%   Xin  - resulting extended matrix
%   V  - right subspace
%   kk - truncation rank for resulting SVD matricies
%   ||A - U*S*V'|| / ||A|| < tol
%
% -------------------------------------------------------------------------

nS = length(idxS);

% find the optimal size for DEIM selection matrix
coeffdeim = min([size(Uin,2), 3*numsv]);
Uin = Uin(:,1:coeffdeim);

% construct the selection matrix Pin for Uin basis
% it is important to note, that DEIM has to be applied 
% prior to 

% check if DEIM++ is available
if (3 == exist('mexDeim','file'))
    [~, Pin] = mexDeim(Uin);
else
    [~, Pin] = deim(Uin);
end


%%
% transform the DEIM points into DEIM voxels: i.e. add the 3 components of
% the selected points

% get S domain coordinates
xs = r(:,:,:,1); xs = xs(idxS);
ys = r(:,:,:,2); ys = ys(idxS);   
zs = r(:,:,:,3); zs = zs(idxS);

% find all the voxels with at least one component selected by deim
Pin3D = [Pin(1:nS,:), Pin(nS+1:2*nS,:), Pin(2*nS+1:3*nS,:)];
vec3D = sum(Pin3D,2);
idxD = find(vec3D);
idxD = sort(idxD);

% get deim points coordinates
xds = xs(idxD);
yds = ys(idxD);
zds = zs(idxD);
Dcoord = [xds(:), yds(:), zds(:)];

% 3D idx and build the Pin for the DEIM approach
idxD = [idxD; nS+idxD; 2*nS+idxD];
e = speye(3*nS);
Pin = sparse(e(:, idxD));
ndeim = size(Pin,2);

% Obtain the DEIM extended inverse matrix
Xin= (Pin.'*Uin)\speye(ndeim,ndeim);