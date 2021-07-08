function [Q] = compute_RWG_moments(RWG, GL_order, i_sie, order) 
% function [Q] = compute_RWG_moments(RWG, order) 
% computes RWG moments of orders <= order
% first version of this function allows to compute 
% orders 0 (3 equations) and 1 (3 + 21 equations);


switch order 
    case 0
        % compute coefficients for 0 order expansion
        Q = current_conservation_RWG(RWG, GL_order,i_sie);
        
    case 1
        % allocate memory for 0,1-st order expansion functions 
        Q = zeros(12,1);
        
        % compute coefficients for order <= 1
        Q(1:3)  = current_conservation_RWG(RWG, GL_order);
        Q(4:12) = first_moment_RWG(RWG, GL_order);
         
    otherwise
                       
        error('current compute_RWG_moments() version computes moments up to 1-st order!\n');
end;


