function [I] = current_conservation_PWX(r_center, res, GL_order, PWX_flag)
% function [I] = current_conservation_PWX(r_center, res, GL_order, PWL_flag)


%% get number of basis components for each coordinate
N_q = pwx_size(0);
N_l = pwx_size(PWX_flag) / pwx_size(0);
 
% get problem dimensions
[m,n] = size(r_center);

% number of voxels along one direction
L = round(n.^(1 / 3));

% Allocte memory 
I       = zeros(m,n * N_l * N_q); 
alpha   = zeros(N_l,1); 

%% get quadrature points

[wt_vie,z] = gauss_1d_sie(GL_order);


%% compute entries and fill the matrix

for basis = 1:N_l


    %integrate along x-axis
    for i1_vie = 1:length(wt_vie)
        % integrate along y-axis
        for j2_vie = 1:length(wt_vie)
            % integrate along z - axis
            for k3_vie = 1:length(wt_vie)
                
                Weight       = wt_vie(i1_vie) * wt_vie(j2_vie) * wt_vie(k3_vie) / 8.0;
                alpha(basis) = alpha(basis) + Weight * pwx_fct(z(i1_vie), z(j2_vie), z(k3_vie), basis);
                
            end;
        end;
    end;
        
end;

idx_0x = 1:n;
idx_0y = (n * N_l + 1):( n * (N_l + 1));
idx_0z = (n * (N_l + 1) + 1):( n * (N_l + 2));

% fill in constant terms
I(1,idx_0x) = ones(1,n) .* alpha(1) * res^3; 
I(2,idx_0y) = ones(1,n) .* alpha(1) * res^3;
I(3,idx_0z) = ones(1,n) .* alpha(1) * res^3;

%fill linear terms if applicable
if(PWX_flag)
    
    % set indicies for linear components for q = x;
    idx_1x_x = (n + 1):(2 * n); 
    idx_1x_y = (2 * n + 1):(3 * n);
    idx_1x_z = (3 * n + 1):(4 * n);
    
    % set indicies for linear components for q = y;
    idx_1y_x = (n * (N_l + 1) + 1):(n * (N_l + 2)); 
    idx_1y_y = (n * (N_l + 2) + 1):(n * (N_l + 3));
    idx_1y_z = (n * (N_l + 3) + 1):(n * (N_l + 4));
    
    % set indicies for linear components for q = z;
    idx_1z_x = (n * (2 * N_l + 1) + 1):(n * (2 * N_l + 2)); 
    idx_1z_y = (n * (2 * N_l + 2) + 1):(n * (2 * N_l + 3));
    idx_1z_z = (n * (2 * N_l + 3) + 1):(n * (2 * N_l + 4));
    
    % fill matrix
    I(1,idx_1x_x) = ones(1,n) * alpha(2) * res^3;
    I(1,idx_1x_y) = ones(1,n) * alpha(3) * res^3;
    I(1,idx_1x_z) = ones(1,n) * alpha(4) * res^3;

    % fill matrix
    I(2,idx_1y_x) = ones(1,n) * alpha(2) * res^3;
    I(2,idx_1y_y) = ones(1,n) * alpha(3) * res^3;
    I(2,idx_1y_z) = ones(1,n) * alpha(4) * res^3;
    
    % fill matrix
    I(3,idx_1z_x) = ones(1,n) * alpha(2) * res^3;
    I(3,idx_1z_y) = ones(1,n) * alpha(3) * res^3;
    I(3,idx_1z_z) = ones(1,n) * alpha(4) * res^3;
end;

% I =   


