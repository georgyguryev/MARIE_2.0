function [rhbm_out] = RHBM_4comp_cuboid_(obj,res,Len,Cnt,e_r,s_e,dens)

    % Generates a four hmogeneous compartment cuboid phantom
    % Author: Ilias I. Giannakopoulos, Cambridge, MA, 2019
    
    rhbm_out.name = '4_Compartment_Cuboid';

    r = obj.RHBM_generatedomain(res,Len(1)/2,Len(2)/2,Len(3)/2);

    [n1,n2,n3,~] = size(r);
    Cnt = squeeze(Cnt);
    Len = squeeze(Len);
    
    epsilon_r = ones(n1,n2,n3);
    sigma_e = zeros(n1,n2,n3);
    rho = zeros(n1,n2,n3);
    
    cube = @(r)( ( abs(r(:,:,:,1)-Cnt(1)) <= Len(1)/2 ) & ( abs(r(:,:,:,2)-Cnt(2)) <= Len(2)/2 ) & ( abs(r(:,:,:,3)-Cnt(3)) <= Len(3)/2 ) ) ;
    pointsI = cube(r); 
    idxS = find(pointsI(:)); 
    
    cube_back_left = @(r)( ( (r(:,:,:,1)-Cnt(1)) <= 0 ) & ((r(:,:,:,1)-Cnt(1)) >= -Len(1)/2 ) & ( (r(:,:,:,2)-Cnt(2)) <= 0 ) & ( (r(:,:,:,2)-Cnt(2)) >= -Len(2)/2 ) & ( abs(r(:,:,:,3)-Cnt(3)) <= Len(3)/2 ) ) ;
    pointsI = cube_back_left(r); 
    idxS1 = find(pointsI(:)); 
    
    cube_front_left = @(r)( ( (r(:,:,:,1)-Cnt(1)) <= 0 ) & ((r(:,:,:,1)-Cnt(1)) >= -Len(1)/2 ) & ( (r(:,:,:,2)-Cnt(2)) <= Len(2)/2 ) & ( (r(:,:,:,2)-Cnt(2)) >= 0 ) & ( abs(r(:,:,:,3)-Cnt(3)) <= Len(3)/2 ) ) ;
    pointsI = cube_front_left(r); 
    idxS2 = find(pointsI(:)); 
    
    cube_back_right = @(r)( ( (r(:,:,:,1)-Cnt(1)) <= Len(1)/2 ) & ((r(:,:,:,1)-Cnt(1)) >= 0 ) & ( (r(:,:,:,2)-Cnt(2)) <= 0 ) & ( (r(:,:,:,2)-Cnt(2)) >= -Len(2)/2 ) & ( abs(r(:,:,:,3)-Cnt(3)) <= Len(3)/2 ) ) ;
    pointsI = cube_back_right(r); 
    idxS3 = find(pointsI(:)); 

    cube_front_right = @(r)( ( (r(:,:,:,1)-Cnt(1)) <= Len(1)/2 ) & ((r(:,:,:,1)-Cnt(1)) >= 0 ) & ( (r(:,:,:,2)-Cnt(2)) <= Len(2)/2 ) & ( (r(:,:,:,2)-Cnt(2)) >= 0 ) & ( abs(r(:,:,:,3)-Cnt(3)) <= Len(3)/2 ) ) ;
    pointsI = cube_front_right(r); 
    idxS4 = find(pointsI(:)); 
    
    epsilon_r(idxS1) = e_r(1); 
    epsilon_r(idxS2) = e_r(2); 
    epsilon_r(idxS3) = e_r(3); 
    epsilon_r(idxS4) = e_r(4); 

    sigma_e(idxS1) = s_e(1);
    sigma_e(idxS2) = s_e(2);
    sigma_e(idxS3) = s_e(3);
    sigma_e(idxS4) = s_e(4);

    rho(idxS1) = dens(1);
    rho(idxS2) = dens(2);
    rho(idxS3) = dens(3);
    rho(idxS4) = dens(4);
    
    rhbm_out.epsilon_r = epsilon_r;
    rhbm_out.sigma_e = sigma_e;
    rhbm_out.rho = rho;
    rhbm_out.r = r;
    rhbm_out.idxS = idxS;
%     rhbm_out.idxS1 = idxS1;
%     rhbm_out.idxS2 = idxS2;
%     rhbm_out.idxS3 = idxS3;
%     rhbm_out.idxS4 = idxS4;
    
end

