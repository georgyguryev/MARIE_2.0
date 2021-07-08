function [rhbm_out] = RHBM_cuboid_(obj,res,Len,Cnt,e_r,s_e,dens)

    % Generates a homogeneous cuboid phantom
    % Author: Ilias I. Giannakopoulos, Cambridge, MA, 2019
    
    rhbm_out.name = 'Homogeneous_Cuboid';

    r = obj.RHBM_generatedomain(res,Len(1)/2,Len(2)/2,Len(3)/2);

    [n1,n2,n3,~] = size(r);
    Cnt = squeeze(Cnt);
    Len = squeeze(Len);
    
    cube = @(r)( ( abs(r(:,:,:,1)-Cnt(1)) <= Len(1)/2 ) & ( abs(r(:,:,:,2)-Cnt(2)) <= Len(2)/2 ) & ( abs(r(:,:,:,3)-Cnt(3)) <= Len(3)/2 ) ) ;
    pointsI= cube(r); 
    idxS = find(pointsI(:)); 
    
    epsilon_r = ones(n1,n2,n3);
    epsilon_r(idxS) = e_r; 
    rhbm_out.epsilon_r = epsilon_r;

    sigma_e = zeros(n1,n2,n3);
    sigma_e(idxS) = s_e;
    rhbm_out.sigma_e = sigma_e;

    rho = zeros(n1,n2,n3);
    rho(idxS) = dens;
    rhbm_out.rho = rho;
    
    rhbm_out.idxS = idxS;
    rhbm_out.r = r;
    
end
