function    [E] =  E_DGF_pwc_Col(r_center, res, r_col, GL_order, f) 
% function  [E] =  E_DGF_pwc_Col(r_center, res, r_col, GL_order, f) 
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


%% get size of output vector
[rows,cols] = size(r_col);

[~, N]  = size(r_center);
M       = rows * cols;

E = zeros(M, 3 * N);


%% EM constants
mu = 4*pi*1e-7;
co = 299792458;
% eo = 1/co^2/mu;

% frequency-dependend parameters
lambda = co / f;
omega  = 2 * pi * f;
ko     = 2 * pi / lambda; 

%% get SIE points 

[wt_vie,z] = gauss_1d_sie(GL_order);

%% this part is developed specificaly for a single PWC sopport

for voxel = 1:N

    for point = 1:cols
        
        % get current collocation point
        r_cur = r_col(:,point);
        
        % create current electric field components, init with zeros
        E_cur = zeros(3);
        
        
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
                    Q     = (ko^2 * R^2  - 1i * ko * R - 1) / P;
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
                    
                    
                    E_cur = E_cur + Weight .* mult .* res^3 .* DGF;
                    
                end;
            end;
        end;
        
        % multiply by scaling factor
        E_cur = 1i .* omega .* 1/(ko^2) .* mu .* E_cur;
        
        E(rows*(point-1) + 1:rows*point, rows*(voxel-1) + 1:rows*voxel) =  E_cur;

    end;

end;
