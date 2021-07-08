function rank = find_max_rank(S, tol)

[~,n] = size(S);

ranks = zeros(n,1);

for i = 1:n
    ranks(i) = find_rank_from_S(S(:,i), tol);
end

rank = max(ranks);

