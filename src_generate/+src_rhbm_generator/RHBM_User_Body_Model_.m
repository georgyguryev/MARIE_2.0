function [rhbm_out] = RHBM_User_Body_Model_(obj)

    % Generates a user defined realistic body model
    % If you wish to load your own body model from the graphical user
    % interface, feel free to modify this function as you want (now it
    % creates an air sphere). 
    % The function should return a struct rhbm_out with fields
    % r, epsilon_r, sigma_e, rho, idxS, name
    % Author: Your name

    rhbm_out.name = 'User_body_model';

    res = 5e-3;
    Rad = 0.15;
    Cnt = [0 0 0];
    e_r = 2;
    s_e = 2;
    dens = 1;
    
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

