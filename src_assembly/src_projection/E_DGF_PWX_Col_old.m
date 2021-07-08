function    [E] =  E_DGF_PWX_Col_old(r_center, res, r_col, GL_order, f, PWL_flag) 
% function  [E] =  E_DGF_PWX_Col(r_center, res, r_col, GL_order, f, PWL_flag) 
% this function computes three components of the electric fields 
% created by expansion voxels at collocation locations r_col
%
% -----------------------------------------------------------------------
%                      CAUTION! 
% -----------------------------------------------------------------------
% The resulting vector E is the actual electric field 
% The field components in vector E are stored in format: 
%
% E = [[E_1x(r_1), E_1y(r1), E_1z(r1), ..., E_nx(r1), E_ny(r1), E_nz(r1)]; 
%      [  ... ];
%     [E_1x(r_m), E_1y(rm), E_1z(rm), ..., E_nx(rm), E_ny(rm), E_nz(rm)]].'         
%
% -----------------------------------------------------------------------
%
%
% -----------------------------------------------------------------------
%                      INPUT:
% -----------------------------------------------------------------------
%       r_center - 
%       res      -
%       r_col    - 
%       GL_order - 
%       f        -
%       PWL_flag - 1 for PWL expansion basis, 0 - PWC exp. basis    

%% get size of output vector

pwx_sz  = pwx_size(PWL_flag);
N_basis = pwx_sz / 3;

% rows - 3 coordinate components; cols - number of collocation voxels
[rows,cols] = size(r_col);

[~, N]      = size(r_center);
M           = rows * cols;

% allocate memory for output matrix
E = zeros(M, pwx_sz * N);


%% EM constants
mu = 4*pi*1e-7;
co = 299792458;
% eo = 1/co^2/mu;

% frequency-dependend parameters
lambda = co / f;
omega  = 2 * pi * f;
ko     = 2 * pi / lambda; 

%% get 1D VIE points 

[wt_vie,z] = gauss_1d_sie(GL_order);

%% this part is developed specificaly for a single PWC sopport

for voxel = 1:N

    for point = 1:cols
        
        % get current collocation point
        r_cur = r_col(:,point);
        
        E_cur = zeros(3, 3 * N_basis);

        for  basis = 1:N_basis
        
            % create current electric field components, init with zeros
            
            
            % Integrate PWC support
            
            % along x axis
            for i1_vie = 1:length(wt_vie)
                
                x_src = r_center(1, voxel) + res / 2 * z(i1_vie);
                
                % along y axis
                for j2_vie = 1:length(wt_vie)
                    
                    y_src = r_center(2, voxel) + res / 2 * z(j2_vie);
                    
                    % along z axis
                    for k3_vie = 1:length(wt_vie)
                        
                        z_src = r_center(3, voxel) + res / 2 * z(k3_vie);
                        
                        % get coordinates of current source
                        r_src = [x_src; y_src; z_src];
                        
                        % define resulting quadrature weight
                        Weight = wt_vie(i1_vie) * wt_vie(j2_vie) * wt_vie(k3_vie) / 8.0;
                        
                        % find vector and distance from current source point to
                        % collocation point
                        R_vec =  r_src - r_cur;
                        R     =  norm(R_vec);
                        
                        % define multiplication coeeficients
                        kappa = exp(-1i * ko * R) / (4 * pi * R^3);
                        P     = (3 + 3 * 1i * ko * R - ko^2 * R^2) / R^2;
                        Q     = (1 + 1i * ko * R - ko^2 * R^2) / P;
                        mult  = kappa * P;
                        
                        % compute DGF entries
                        Gxx   = R_vec(1) * R_vec(1) - Q;
                        Gxy   = R_vec(1) * R_vec(2);
                        Gxz   = R_vec(1) * R_vec(3);
                        Gyy   = R_vec(2) * R_vec(2) - Q;
                        Gyz   = R_vec(2) * R_vec(3);
                        Gzz   = R_vec(3) * R_vec(3) - Q;
                        
                        % form dyadic Green's function
                        DGF = [Gxx, Gxy, Gxz;
                               Gxy, Gyy, Gyz;
                               Gxz, Gyz, Gzz];
                        
                        
                        E_cur(:,3*(basis-1)+1:3*basis) = E_cur(:,3*(basis-1)+1:3*basis) + ...
                                                         Weight .* mult .* res^3 .* DGF * ...
                                                         pwx_fct(z(i1_vie),z(j2_vie),z(k3_vie), basis);
                        
                    end
                end
            end
            
        end
        
        % multiply by scaling factor
        E_cur = 1i .* omega .* 1/(ko^2) .* mu .* E_cur;
        
        
        E(rows*(point-1) + 1:rows*point, N_basis*rows*(voxel-1)+1:N_basis*rows*voxel) = E_cur;

    end;

end;
