function [rhbm_out] = RHBM_sphere_(obj,res,Rad,Cnt,e_r,s_e,dens)

    % Generates a homogeneous spherical phantom
    % Author: Ilias I. Giannakopoulos, Cambridge, MA, 2019
    
    rhbm_out.name = 'Homogeneous_Sphere';

    r = obj.RHBM_generatedomain(res,Rad,Rad,Rad);

    [n1,n2,n3,~] = size(r);
    Cnt = squeeze(Cnt);
    
    sphere = @(r)( (r(:,:,:,1) - Cnt(1) ).^2 + ( r(:,:,:,2) - Cnt(2) ).^2 + ( r(:,:,:,3) - Cnt(3) ).^2 < Rad^2) ;
    pointsI= sphere(r); 
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

