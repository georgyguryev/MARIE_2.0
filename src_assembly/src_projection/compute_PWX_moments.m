function [Q] = compute_PWX_moments(vox_centers,RWG, res, GL_order, order, PWX_flag) 
% function [Q] = compute_PWC_moments(RWG, order) 
% computes PWC moments of orders <= order
% first version of this function allows to compute 
% orders 0 (3 equations) and 1 (3 + 21 equations);

[m,n] = size(vox_centers);


switch order 
    case 0
        % compute coefficients for 0 order expansion
        Q = current_conservation_PWX(vox_centers, res, GL_order, PWX_flag);
        
    case 1
        
        if(1 == PWX_flag)
            N_basis = 4;
        else
            N_basis = 1;
        end;
            
        
        % allocate memory for 0,1-st order expansion functions 
        Q = zeros(12, m * n * N_basis);
        
        % compute coefficients for order <= 1
        Q(1:3,:)  = current_conservation_PWX(vox_centers, res, GL_order, PWX_flag);
        Q(4:12,:) = first_moment_PWX(vox_centers, RWG, res, GL_order, PWX_flag);
        
         
    otherwise
                       
        error('current compute_RWG_moments() version computes moments up to 1-st order!\n');
end;


