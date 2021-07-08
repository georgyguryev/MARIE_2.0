function [U,S,V,kk] = rSVD_MRGF(Fdirect,Fadjoint,r,idxJ,idxE,Lmax,tol,blocksize)
%% Computes the (randomized) SVD of an operator

%   INPUT:
%   Fdirect  - handle to direct operator
%   Fadjoint - handle to the adjoint operator
%   r - 4D tensor of voxel coordinates, used to extract [L,M,N] 
%   m - Left dimension of Fdirect
%   n - right dimension of Fdirect
%   idxJ - 1D index of voxels used for dipole excitation
%   idxE - 1D index of voxels used for 
%   Lmax - maximum number of random excitations (if <1 automatic selection)
%   tol - tolerance for truncation
%   blocksize - block for iterative approaches
%
%   OUTPUT:
%   U  - left subspace
%   S  - singular values
%   V  - right subspace
%   kk - truncation rank for resulting SVD matricies
%   ||A - U*S*V'|| / ||A|| < tol
%
% -------------------------------------------------------------------------

if (nargin < 5)
    fprintf(1, '/nERROR: invalid arguments\n');
    U = [];
    S = [];
    V = [];
    return
end

if isempty(Lmax) || (nargin < 6)
    Lmax = min(500,N);
end

if isempty(tol) || (nargin < 7)
    tol = 1e-3;
end

if isempty(blocksize) || (nargin < 8)
    blocksize = 100;
end


% -------------------------------------------------------------------------
%   Random excitations
%   Gm: (mxl)
%   Gn: (nxl)
%
%   Action on random matrices
%   AGn = A*Gn: (mxl) 
%   AGm = A'*Gm: (nxl)
%
%  ON basis
%  [Q,~] = qr(AGn);    
%  [W,~] = qr(AGm);
%
%
%  T = Q'*A*W;
% -------------------------------------------------------------------------

m = size(idxE,1);
n = size(idxJ,1);

% -------------------------------------------------------------------------
%            Generate the random excitations
% -------------------------------------------------------------------------

fid = 1;
tini = tic;

% if on-the-fly control of SV drop
if (Lmax < 1)
    % generate the random excitation matrix, and initialize the fields
    fprintf(fid, '\n\n RSVD_USV starting with on-the-fly control and tol = %1.1d',tol);
else
    % generate the random excitation matrix, and initialize the fields
    fprintf(fid, '\n\n RSVD_USV starting with %d excitations', Lmax);
end


% -------------------------------------------------------------------------
%            Q and W
% -------------------------------------------------------------------------
fprintf(fid, '\n\n-------------------------- \n');
fprintf(fid, 'ON basis (Q) for left subspace \n');

t2 = tic;
[Q,Nrand_Q] = getQ_MRGF(Fdirect,r,m,n,idxJ,idxE,Lmax,tol,blocksize);
timeQR_left = toc(t2);
fprintf(fid,'\n\nelapsed time %g, with %d excitations', timeQR_left,Nrand_Q);

fprintf(fid, '\n\n-------------------------- \n');
fprintf(fid, 'ON basis (W) for right subspace \n');

t2=tic;
[W,Nrand_W] = getQ_MRGF(Fadjoint,r,n,m,idxE,idxJ,Lmax,tol,blocksize);
timeQR_right = toc(t2);
fprintf(fid,'\n\nelapsed time %g, with %d excitations', timeQR_right,Nrand_W);
%%

[L,M,N,~] = size(r);
% -------------------------------------------------------------------------
%  T = Q'*A*W;
% -------------------------------------------------------------------------
t1 = tic;

% check minimum columns between Q and W

if size(Q,2) < size(W,2)
    %  compute  A'*Q
    U = zeros(n,size(Q,2));
    for i = 1:size(Q,2)
        
        % form excitation tensor to apply operator
        ve_LMN3 = zeros(L,M,N,3);
       
        % form current excitation vector 
        ve = Q(:,i); 
        
        % map current excitation vector to [L,M,N,3] format
        ve_LMN3(idxE) = ve;
        
        % adjoint operator
        U_LMN3 = Fadjoint(ve_LMN3);
        U(:,i) = U_LMN3(idxJ);
    end
    % compute T = transpose(A*conj(Q)) * W
    T = U.'*W;
else
    %  compute  A*W
    U = zeros(m,size(W,2));
    for i = 1:size(W,2)
        
        % form excitation tensor to apply operator
        ve_LMN3 = zeros(L,M,N,3);
        
        % form current excitation vector
        ve = W(:,i); 
        
        % map current excitation vector to [L,M,N,3] format
        ve_LMN3(idxJ) = ve;
        
        % adjoint operator
        U_LMN3 = Fdirect(ve_LMN3);
        U(:,i) = U_LMN3(idxE);
    end
    % compute T = ctranspose(Q) * (A*W);
    T = Q'*U;
end
timeinc = toc(t1);
fprintf(fid, '\n\n-------------------------- \n');
fprintf(fid,'T computed, elapsed time %g \n', timeinc);

% -------------------------------------------------------------------------
%          [Ur,S,Vr] = svd(T);
% -------------------------------------------------------------------------

t1 = tic;

% -------------------------------------------------------------------------
% apply SVD and truncation with defined tolerance
% -------------------------------------------------------------------------
[U, S, V] = svd(T,'econ');

for kk = 2:length(S)
    if S(kk,kk) < tol*S(1,1)
        Nsv = kk;
        break
    end
end

if kk == length(S)
    Nsv = kk;
end

% -------------------------------------------------------------------------
% Final output
% -------------------------------------------------------------------------
U = Q*U; 
V = W*V;

timeinc = toc(t1);
fprintf(fid,'\n SVD(T) done, time %g', timeinc);

fprintf(fid, '\n\n-------------------------- \n');
fprintf(fid, 'Random SVD operation done\n');
fprintf(fid, ' elapsed time %g seconds\n',toc(tini));
fprintf(fid, ' #Singular values  %d\n',Nsv);

end