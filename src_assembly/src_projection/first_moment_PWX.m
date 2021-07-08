function [I] = first_moment_PWX(r_center, RWG, res, GL_order, PWL_flag)

% number of linear moments (for x,y and z) 
moments = 3;

%% get number of basis components for each coordinate

if(1 == PWL_flag)
    N_basis = 4;        % number of basis vector components
else
    N_basis = 1;        % number of basis vector components
end;

% get problem dimensions
[m,n] = size(r_center);

% Allocte memory 
I     = zeros(m * moments, m * n * N_basis);

%% get center of RWG common edge

Ct = RWG.Ct;

%% get quadrature points

[wt_vie,z] = gauss_1d_sie(GL_order);


%% compute entries and fill the matrix

% loop over voxels
for vox = 1:n
    
    %get current voxel's center
    r_c = r_center(:,vox);
    
    % loop over first moment entries
    for moment = 1:moments
        
        % allocate buffer to store coefficients
        sub_mat = zeros(m, m * N_basis);
        buf     = zeros(N_basis,1);
        
        
        buf(1)  = res^3 * (r_c(moment) - Ct(moment));
        
        % 
        for basis = 2:N_basis
            % nonzero elements occur only when moment type matches basis
            % type (i.e. moment for x component and basis el. for x comp.)
            if((basis - 1) == moment)
                
                % compute matching component with GL quadrature
                for t_vie = 1:length(wt_vie)
                    
                    x = res / 2 * z(t_vie) + r_c(basis-1);
                    buf(basis) = buf(basis) + wt_vie(t_vie) * (x - Ct(moment)) * (x - r_c(basis-1));
                end;
                
                scaling_fct = res^2 / 2;                
                buf(basis) = buf(basis) * scaling_fct;

            end;
        end;
        
        sub_mat = kron(buf, eye(m)).';
        
        row_idx = m * (moment - 1)+1:m * moment;
        col_idx = m * N_basis * (vox -1) + 1: m * N_basis * vox;
        
        I(row_idx, col_idx) = sub_mat;
    end;
    
end;



