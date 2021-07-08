function [Q] = compute_PWC_moments(vox_centers,RWG, res, GL_order, order) 
% function [Q] = compute_PWC_moments(RWG, order) 
% computes PWC moments of orders <= order
% first version of this function allows to compute 
% orders 0 (3 equations) and 1 (3 + 21 equations);

[m,n] = size(vox_centers);


switch order 
    case 0
        % compute coefficients for 0 order expansion
        Q = current_conservation_PWC(vox_centers, res, GL_order);
        
    case 1
        % allocate memory for 0,1-st order expansion functions 
        Q = zeros(12, m * n);
        
        % compute coefficients for order <= 1
        Q(1:3,:)  = current_conservation_PWC(vox_centers, res, GL_order);
        Q(4:12,:) = first_moment_PWC(vox_centers, RWG, res, GL_order);
        
         
    otherwise
                       
        error('current compute_RWG_moments() version computes moments up to 1-st order!\n');
end;


