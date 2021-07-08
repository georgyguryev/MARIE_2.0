function [rank] = find_rank_from_S(s, tol)
% function find_rank_cross_2d() finds rank of cross approximation based on
% relative  of a submatrix of maximal volume

    % total number of singular values
    N_s = length(s);
    
    s_2 = s.^2;
    
    % compute a square of Frobenius norm
%     full_squares_sum = sum(s_2);
    % initialize partial square sum
    partial_square_sum = 0;
    
    % find current rank
    for i = N_s:-1:1
        % compute relative truncation error (Frobenius norm)
        partial_square_sum = partial_square_sum + s_2(i);
        rel_trunc_error = sqrt(partial_square_sum) / s(1);
        % stop when truncation error exceeds the tolerance
        if (rel_trunc_error > tol)
            break;
        end
    end
    rank = i;
end