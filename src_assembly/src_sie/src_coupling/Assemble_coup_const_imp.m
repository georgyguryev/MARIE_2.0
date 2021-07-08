function [Zbc_const] = Assemble_coup_const_imp(R1,R2,R3,RO,NO,IDX,L1,ko,res,Np_2D,Z1,Z2,Z3,wp,GL_VIE_order, idx_imp_near)
% function Assemble_coup_lin(Scoord,index,etod,node,elem,freq,LEVEL_DVrule,i_sie)
% assebles coupling matrix between descritized coil surface (RWG) and a
% tissue (with piece-wise constant field approximation, PWL); function is
% designed to construct a coupling vector(matrix) for given dof (dofs)
%
% -----------------------------------------------------------------------
% INPUT:  
%         R1  - first vertices of parent triangles 
%         R2  - second vertices of parent triangles
%         R3  - third vertices of parent triangles
%         RO  - 3D coordinates of body voxel centers
%         NO  - number of bosy voxels
%         IDX - 
%         L1  - length of an inner edge
%         ko  - free space wavenumber
%         res - resolution of the body model
%         Np_2d - number of Dunavant quadrature points 
%         Z1, Z2, Z3 - coordinates of dunavant quadrature points
%         wp  - weights for Dunavant quadrature rule
%         GL_VIE_order - order of Gauss-Legendre 1D quadratures
%       
% ----------------------------------------------------------------------- 
% OUTPUT:
%       - Zbc_lin - coupling matrix with between RWG and PWL bases
% ----------------------------------------------------------------------- 

%% get 1D quadrature points for VIE
[wt_vie_far,z_vie_far]       = gauss_1d_sie(GL_VIE_order.far);
% [wt_vie_medium,z_vie_medium] = gauss_1d_sie(GL_VIE_order.medium);
[wt_vie_near,z_vie_near]     = gauss_1d_sie(GL_VIE_order.near);


%% get sizes 
Nq   = pwx_size(0);
Nb   = NO * Nq;
Ns   = size(R1,2);        % number of supporting triangles

%% allocate output matrix

Zbc_const = zeros(NO, Nq);

%%
RR = [R1, R2, R3];

%% assemble 

for vox = 1:NO
    
    
    if isempty(find(idx_imp_near == vox))
        wt_vie = wt_vie_far;
        z_vie  = z_vie_far;
    else
        wt_vie = wt_vie_near;
        z_vie  = z_vie_near;
    end
    
    Nvie  = length(z_vie);

    
    Integral_1 = zeros(Nq,1);
    
    RO_c = RO(:,vox);
    
    for tr = 1:Ns
        
        % find index of current opposite vertice
        vert_idx = find(IDX(:,tr));
        
        % get coordinates of the opposite node
        r1 = RR(:,(vert_idx - 1) * 2  + tr);
        
        % integrate over the surface
        for i_sie = 1:Np_2D
            
            % source vector (on a coil surface)
            r_src  = Z1(i_sie) * R1(:,tr) + Z2(i_sie) * R2(:,tr) + Z3(i_sie) * R3(:,tr);
            wp_sie = wp(i_sie);
            
            % form observation coordinate x_vie
            for i_vie = 1:Nvie
                %                        r_obs(1) =  RO(1,vox) +  z_vie(i_vie) * res / 2;
                x        = z_vie(i_vie);
                wt_vie_x = wt_vie(i_vie);
                
                % form observation coordinate y_vie
                for j_vie = 1:Nvie
                    %                            r_obs(2) = RO(2,vox) + z_vie(j_vie) * res / 2;
                    y        = z_vie(j_vie);
                    wt_vie_y = wt_vie(j_vie);
                    
                    % form observation coordinate z_vie
                    for k_vie = 1:Nvie
                        
                        z        =  z_vie(k_vie);
                        wt_vie_z = wt_vie(k_vie);
                        
                        r_obs = RO_c + [x; y; z] .* res ./ 2;
                        %                                r_obs(3) = RO(3,vox) + z_vie(k_vie) * res / 2;
                        
                        % form distance vector R
                        R_vec = r_obs - r_src;
                        
                        % form vector rho
                        rho = (r_src - r1) * (1.0 - 2.0 * (tr - 1));
                        
                        % define resulting quadrature weight
                        Weight = wt_vie_x * wt_vie_y * wt_vie_z * wp_sie / 8.0;
                        
                        % distance
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
                                                
                        
                        Integral_1 = Integral_1 + Weight .* mult .* DGF * rho;
                        
                    end
                end
            end
        end
    end
    
    Zbc_const(vox,:) = Integral_1 * L1;

end


%% reshape Coupling vector

Zbc_const = reshape(Zbc_const, Nb, 1);
