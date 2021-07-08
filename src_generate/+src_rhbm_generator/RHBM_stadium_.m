function [rhbm_out] = RHBM_stadium_(obj,res,Rad,side,Len,Cnt,e_r,s_e,dens)

    % Generates a homogeneous stadium phantom
    % Author: Ilias I. Giannakopoulos, Cambridge, MA, 2019

    rhbm_out.name = 'Homogeneous_Stadium';

    r = obj.RHBM_generatedomain(res,Rad,Rad+side/2,Len/2);

    [n1,n2,n3,~] = size(r);
    Cnt = squeeze(Cnt);
    
    Cnt_1  = Cnt + [0 side/2 0];
    Cnt_2  = Cnt - [0 side/2 0];
    
    semi_cylinder_1 = @(r)( (abs(r(:,:,:,2)-Cnt_1(2) + 1i*(r(:,:,:,1)-Cnt_1(1)) ) <= Rad) & ( r(:,:,:,2)-Cnt_1(2) >= 0 ) & ( (r(:,:,:,3)-Cnt_1(3)) <= ( Len/2 ) ) & ( (r(:,:,:,3)-Cnt_1(3)) >= (-Len/2)) );
    points_1 = semi_cylinder_1(r); 
    idxS_1 = find(points_1(:)); 
    
    semi_cylinder_2 = @(r)( ( abs(r(:,:,:,2)-Cnt_2(2) + 1i*(r(:,:,:,1)-Cnt_2(1)) ) <= Rad) & ( r(:,:,:,2)-Cnt_1(2) <= 0 ) & ( (r(:,:,:,3)-Cnt_2(3)) <= ( Len/2 ) ) & ( (r(:,:,:,3)-Cnt_2(3)) >= (-Len/2)) );
    points_2 = semi_cylinder_2(r); 
    idxS_2 = find(points_2(:)); 
    
    cuboid = @(r)( ( abs(r(:,:,:,1)-Cnt(1)) <= Rad ) & ( abs(r(:,:,:,2)-Cnt(2)) <= side/2 ) & ( abs(r(:,:,:,3)-Cnt(3)) <= Len/2 ) ) ;
    points_3 = cuboid(r); 
    idxS_3 = find(points_3(:)); 
    
    epsilon_r = ones(n1,n2,n3);
    epsilon_r(idxS_1) = e_r; 
    epsilon_r(idxS_2) = e_r; 
    epsilon_r(idxS_3) = e_r; 
    rhbm_out.epsilon_r = epsilon_r;

    sigma_e = zeros(n1,n2,n3);
    sigma_e(idxS_1) = s_e;
    sigma_e(idxS_2) = s_e;
    sigma_e(idxS_3) = s_e;
    rhbm_out.sigma_e = sigma_e;

    rho = zeros(n1,n2,n3);
    rho(idxS_1) = dens;
    rho(idxS_2) = dens;
    rho(idxS_3) = dens;
    rhbm_out.rho = rho;
    
    rhbm_out.idxS = find(abs(sigma_e(:)));
    rhbm_out.r = r;
    
end

