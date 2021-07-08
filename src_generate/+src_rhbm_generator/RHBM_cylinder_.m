function [rhbm_out] = RHBM_cylinder_(obj,res,Rad,Len,Cnt,e_r,s_e,dens)

    % Generates a homogeneous cylindrical phantom
    % Author: Ilias I. Giannakopoulos, Cambridge, MA, 2019
    
    rhbm_out.name = 'Homogeneous_Cylinder';

    r = obj.RHBM_generatedomain(res,Rad,Rad,Len/2);

    [n1,n2,n3,~] = size(r);
    Cnt = squeeze(Cnt);
    
    cylinder = @(r)( (abs( r(:,:,:,2)-Cnt(2) + 1i*(r(:,:,:,1)-Cnt(1)) ) <= Rad) & ( (r(:,:,:,3)-Cnt(3)) <= ( Len/2 ) ) & ( (r(:,:,:,3)-Cnt(3)) >= (-Len/2)) );
    pointsI= cylinder(r); 
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

