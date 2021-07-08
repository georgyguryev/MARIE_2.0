function P = maxvol(A, tol, N_perm)
% function P = maxvol(A, tol)
% computes indeces of rows that correspond to a submatrix 
% of approximately maximal volume 

% get dimensions
[m,r] = size(A);

% if m < r 
if(m < r)
    [~,r] = size(A');
end

% initialize permutation counter
i_perm = 1;

% compute lu
[L,U,P] = lu(A, 'vector');

% consider only first r indices 
% P = P(1:r);

% get perturbed submatrix A_hat
A_hat = L * U;
A_hat = A_hat(1:r, 1:r);

% compute Up = A * (A_hat)^-1; Up(P,:) should be  == eye(r)
Up = A / A_hat;

Up = Up(P,:);


% find first initial element
max_element = max(max(Up));

while (abs(max_element) > 1 + tol) && (i_perm <= N_perm)
    [i,j] = find(Up == max_element);
    u = Up(:,j);
    u(j) = u(j) - 1;
    u(i) = u(i) + 1;
%     v = Up(P(j),:) - Up(i,:);
    v = Up(j,:) - Up(i,:);
    u = u / Up(i,j);
    Up = Up +  u * v;
    max_element = max(max(Up));
    i_perm = i_perm + 1;
    
    temp = P(j);
    P(j) = P(i);
    P(i) = temp;
end

P = P(1:r).';
