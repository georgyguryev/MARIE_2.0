% Preconditioned GMRES internal iteration routine
%
% Generates Arnoldi relation A V(:,1:m) = V(:,1:m+1) H
% Applies relation (I - C*C') V(:,1:m) = V(:,1:m+1) H if deflation exists
%
% INPUT:  A      N-by-N matrix
%         r      N-by-1 preconditioned residual vector
%         m      number of GMRES iterations to perform
%         L1     first left preconditioner for A
%         L2     second left preconditioner for A
%         R1     first right preconditioner for A
%         R2     second right preconditioner for A
%         C      matrix containing the leading subspace of previous restart
%         tol    specifies the tolerance of the method
% OUTPUT: V      matrix containing orthogonal basis for Krylov subspace 
%         H      upper Hessenburg reduction of matrix operator
%         B      the matrix C'*A*V(:,1:k)
%         k      number of GMRES iterations actually performed
%         resvec vector containing norm of residual at each iteration of GMRES
%
%

function [V,H,B,k,resvec, x] = my_pgmresDR_iter(A,r,m,L1,pL_L,L2,R1,R2,C,tol)

if(isempty(L1))
   existL1 = 0;
else
   existL1 = 1;
  [L1type,L1fun,L1fcnstr] = src_numeric.iterchk(L1);
end
if(isempty(L2))
   existL2 = 0;
else
   existL2 = 1;
   [L2type,L2fun,L2fcnstr] = src_numeric.iterchk(L2);
end
if(isempty(R1))
   existR1 = 0;
else
   existR1 = 1;
  [R1type,R1fun,R1fcnstr] = src_numeric.iterchk(R1);
end
if(isempty(R2))
   existR2 = 0;
else
   existR2 = 1;
   [R2type,R2fun,R2fcnstr] = src_numeric.iterchk(R2);
end

% Determine whether A is a matrix or a function.
[atype,afun,afcnstr] = src_numeric.iterchk(A);

% Preallocate and Initialize V
V = zeros(size(r,1),m+1);
V(:,1) = r / norm(r);

% Preallocate B, H and resvec
B = zeros(size(C,2),m);
H = zeros(m+1,m);
resvec = zeros(m,1);

x = zeros(size(r,1),m+1);


for k = 1:m
    
    w = V(:,k);
    
   % Find w using right preconditioning if available.
   % w = R1 \ ( R2 \ V(:,k) );
   if(existR2) % second right preconditioning
       if strcmp(R2type,'matrix')
           w = R2 \ w;
       else
           w = src_numeric.iterapp('mtimes',R2fun,R2type,R2fcnstr,w);
       end
   end
   if(existR1) % first right preconditioning
       if strcmp(R1type,'matrix')
           
%            w = R1 \ w;
           w = src_numeric.apply_preconditioner(R1, w, []);
       else
           w = src_numeric.iterapp('mtimes',R1fun,R1type,R1fcnstr,w);
       end
   end
   
   % Apply the system operator
   % w = A*w;
   w = src_numeric.iterapp('mtimes',afun,atype,afcnstr,w);
   
   % Apply the Left preconditioning if available.
   % w = L2 \ ( L1 \ w );
   if(existL1) % first left preconditioning
       if strcmp(L1type,'matrix')
%            w = L1 \ w(pL_L);
%            w = L1 \ w;
%            w = L1 .* w;
           
           w = src_numeric.apply_preconditioner(L1, w, pL_L);
       else
           w = src_numeric.iterapp('mtimes',L1fun,L1type,L1fcnstr,w);
       end
   end
   if(existL2) % second left preconditioning
       if strcmp(L2type,'matrix')
%            w = L2 \ w;
%            w = L2 .* w;
           w = src_numeric.apply_preconditioner(L2, w, pL_L);
       else
           w = src_numeric.iterapp('mtimes',L2fun,L2type,L2fcnstr,w);
       end
   end
   
          
   if ~isempty(C)
       % Apply (I-C*C') operator to A*w
       B(:,k) = C' * w;
       w = w - C * B(:,k);
   end
   
   
   % Create next column of V and H
   H(1:k,k) = V(:,1:k)' * w;
   w = w - V(:,1:k)*H(1:k,k);

   H(k+1,k) = norm(w);
   V(:,k+1) = w / H(k+1,k);

   % Initialize right hand side of least-squares system
   rhs = zeros(k+1,1);
   rhs(1) = norm(r);

   % Solve least squares system; Calculate residual norm
   y = H(1:k+1,1:k) \ rhs;
   res = rhs - H(1:k+1,1:k) * y;
   resvec(k) = norm(res);
   
   x(:,k) = V(:,1:k) * y;
   
%    if isstruct(tol)
%        
%        
%        resvec_sie = norm(res(1:dims.N_sie));
%        resvec_vie = norm(res(dims.N_sie + 1:end));
%        
%        if (resvec(k) < tol) && (resvec_sie(k) < tol.sie) &&...
%           (resvec_vie(k) < tol.vie)
%        
%            % truncate preallocated matrices to current size
%            H = H(1:k+1,1:k);
%            V = V(:,1:k+1);
%            B = B(:,1:k);
%            resvec = resvec(1:k);
%            return
%        end
%    end
   
   % check for early convergence
   if resvec(k) < tol
       
       % truncate preallocated matrices to current size
       H = H(1:k+1,1:k);
       V = V(:,1:k+1);
       B = B(:,1:k);
       resvec = resvec(1:k);
       
       return
       
   end
   
end

