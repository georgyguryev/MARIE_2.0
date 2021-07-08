function [phi, P] = deim(U)

% get sizes of the input orthonormal matrix
[n,m] = size(U);

% phi is the vector with the indexes
phi = zeros(m, 1);

% first index: maximum value
[~, phi(1)] = max(abs(U(:, 1)));

% loop on the columns of U
for j = 2:m
    
    idx = phi(1:j-1);
    uL = U(:,j); % select column vector
    c = uL(idx); % select entries of current column vector
    c = U(idx,1:j-1)\c; % compute coeff
    r = uL - U(:,1:j-1)*c; % orthogonalize
    [~, phi(j)] = max(abs(r)); % find maximum entry in column
    
end
% prepare elements for return
phi = sort(phi);
e = speye(n);
P = sparse(e(:, phi));