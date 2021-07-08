function [Q] = first_moment_PWC(vox_centers, RWG, res, GL_order)
% function [Q] = first_moment_PWC(vox_centers, res, GL_order)
% function computes first order PWC moments for 



%% get problem dimensions
[m,n] = size(vox_centers);


% allocate memory for resulting weights
Q = zeros(9, m * n);

%% 
[wt_vie, z] = gauss_1d_sie(GL_order);


ct   = RWG.Ct;


%% assemble matrix of first moments for PWC basis support


for j = 1:n 
    
    % get coordinates of current voxel
    r_c = vox_centers(:,j);
    
    % loop over x, y or z moment
    for i = 1:3
        
        ct_cur  = ct(i);
        r_c_cur = r_c(i);
        
        % scaling factor
        alpha = 0;
        
% this part was commented out since integration of x over -res/2; res/2 is 0        
%         for k_vie = 1:length(wt_vie)
%             
%             alpha = alpha  + res * wt_vie(k_vie) * z(k_vie) / 4;
%             
%         end;
    
        alpha = alpha + r_c_cur - ct_cur;
        
        alpha = alpha * res^3 * eye(3);
        
        Q (3*(i-1)+1:3*i, 3*(j-1)+1:3*j) = alpha;
        
    end;
    
end;