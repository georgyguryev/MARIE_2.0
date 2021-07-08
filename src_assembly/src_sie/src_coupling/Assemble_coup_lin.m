function [Zbc_lin] = Assemble_coup_lin(R1,R2,R3,RO,NO,IDX,L1,ko,res,Np_2D,Z1,Z2,Z3,wp,GL_VIE_order)
% function Assemble_coup_lin(Scoord,index,etod,node,elem,freq,LEVEL_DVrule,i_sie)
% assebles coupling matrix between descritized coil surface (RWG) and a
% tissue (with piece-wise linear field approximation, PWL); function is
% designed to construct a coupling vector(matrix) for given dof (dofs)
%
% -----------------------------------------------------------------------
% INPUT:                                                        
%       
% ----------------------------------------------------------------------- 
% OUTPUT:
%       - Zbc_lin - coupling matrix with between RWG and PWL bases
% ----------------------------------------------------------------------- 

%% get 1D quadrature points for VIE
[wt_vie,z_vie] = gauss_1d_sie(GL_VIE_order);


%% get sizes 
Nb = NO * pwx_size(1);
Nq = pwx_size(0);
Nl = pwx_size(1) / Nq;
Ns = size(R1,2);        % number of supporting triangles
Nvie = length(z_vie);   % get number of quadrature points in VIE 1D 

%% allocate output matrix

Zbc_lin = zeros(NO, Nl, Nq);

%%
RR = [R1, R2, R3];

%% assemble 
for q = 1:Nq
    for l = 1:Nl
        for vox = 1:NO
            
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
                               
                               r_obs = RO_c + [x; y; z] * res / 2;
%                                r_obs(3) = RO(3,vox) + z_vie(k_vie) * res / 2;
                               
                               % form distance vector R
                               R_vec = r_obs - r_src;
                               
                               % form vector rho
                               rho = (r_src - r1) * (1.0 - 2.0 * (tr - 1));
                               
                               % define resulting quadrature weight
                               Weight = wt_vie_x * wt_vie_y * wt_vie_z * wp_sie/ 8.0;
                               
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
                                  
                               pwx_factor = pwx_fct(x,y,z,l);
                        
                               
                               Integral_1 = Integral_1 + Weight .* mult .* pwx_factor .* DGF * rho;
                               
                           end
                       end
                    end                    
                end
            end
            
            Zbc_lin(vox,l,q) = Integral_1(q) * L1;
        end
    end
end

%% reshape Coupling vector

Zbc_lin = reshape(Zbc_lin, Nb, 1);
