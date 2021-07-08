function [rhbm_out] = RHBM_ellipsoid_(obj,res,Rad,Cnt,e_r,s_e,dens)

    % Generates a homogeneous ellipsoid phantom
    % Author: Ilias I. Giannakopoulos, Cambridge, MA, 2019
    
    rhbm_out.name = 'Homogeneous_Ellipsoid';

    r = obj.RHBM_generatedomain(res,Rad(1),Rad(2),Rad(3));

    [n1,n2,n3,~] = size(r);
    Cnt = squeeze(Cnt);
    Rad = squeeze(Rad);
    
    ellipsoid = @(r)( ( (r(:,:,:,1) - Cnt(1))/Rad(1) ).^2 + ( (r(:,:,:,2) - Cnt(2))/Rad(2) ).^2 + ( (r(:,:,:,3) - Cnt(3))/Rad(3) ).^2 < 1) ;
    pointsI= ellipsoid(r); 
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

