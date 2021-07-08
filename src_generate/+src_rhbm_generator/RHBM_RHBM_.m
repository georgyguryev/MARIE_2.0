function [rhbm_out] = RHBM_RHBM_(s,slice)
    
    % Crops a Virtual family model and generates a realistic human body model
    % Author: Ilias I. Giannakopoulos, Cambridge, MA, 2019
    
    rhbm_out.name = strcat(s(1:end-4),'_',slice);
    load(strcat('./src/data/bodies/Virtual_Family/',s),'epsilon_r','sigma_e','r');
    
    %% Cut the air around the whole body
    xd = r(:,:,:,1);
    yd = r(:,:,:,2);
    zd = r(:,:,:,3);
    idxS = find(abs(sigma_e(:)));
    x1  = min(xd(idxS));
    x2 = max(xd(idxS));
    y1  = min(yd(idxS));
    y2 = max(yd(idxS));
    z1  = min(zd(idxS));
    z2 = max(zd(idxS));
    [~,ix1] = min(abs(r(:,1,1,1)-x1));
    [~,ix2] = min(abs(r(:,1,1,1)-x2));
    [~,iy1] = min(abs(r(1,:,1,2)-y1));
    [~,iy2] = min(abs(r(1,:,1,2)-y2));
    [~,iz1] = min(abs(r(1,1,:,3)-z1));
    [~,iz2] = min(abs(r(1,1,:,3)-z2));
    
    %% Reduce the body
    r = r(ix1:ix2,iy1:iy2,iz1:iz2,:);
    epsilon_r = epsilon_r(ix1:ix2,iy1:iy2,iz1:iz2);
    sigma_e = sigma_e(ix1:ix2,iy1:iy2,iz1:iz2);
    idxS = find(abs(sigma_e(:)));

    %% Pick a specific slice
    if strcmp(slice,'Full Body')
        rhbm_out.r = r;
        rhbm_out.epsilon_r = epsilon_r;
        rhbm_out.sigma_e = sigma_e;
        rho = zeros(size(rhbm_out.epsilon_r));
    	rhbm_out.rho = rho;
    	rhbm_out.idxS = idxS;
    elseif strcmp(slice,'Head')
        if strcmp(s(1),'b')
            z2  = -0.52;
        elseif strcmp(s(1),'e')
            z2  = -0.58;
        elseif strcmp(s(1),'d')
            z2  = -0.67;
        end
        [~,iz2] = min(abs(r(1,1,:,3)-z2));
        r= r(:,:,1:iz2,:);
        sigma_e = sigma_e(:,:,1:iz2);
        epsilon_r = epsilon_r(:,:,1:iz2);
    elseif strcmp(slice,'Head and Shoulders')
        if strcmp(s(1),'b')
            z2  = -0.40;
        elseif strcmp(s(1),'e')
            z2  = -0.46;
        elseif strcmp(s(1),'d')
            z2  = -0.53;
        end
        [~,iz2] = min(abs(r(1,1,:,3)-z2));
        r= r(:,:,1:iz2,:);
        sigma_e = sigma_e(:,:,1:iz2);
        epsilon_r = epsilon_r(:,:,1:iz2);
    elseif strcmp(slice,'Right Leg')
        if strcmp(s(1),'b')
            z1 = 0.15;
            x2 = 0;
        elseif strcmp(s(1),'e')
            z1  = 0.1;
            x2 = 0;
        elseif strcmp(s(1),'d')
            z1  = 0.15;
            x2 = 0;
        end
        [~,iz1] = min(abs(r(1,1,:,3)-z1));
        [~,ix2] = min(abs(r(:,1,1,1)-x2));
        r= r(1:ix2,:,iz1:end,:);
        sigma_e = sigma_e(1:ix2,:,iz1:end,:);
        epsilon_r = epsilon_r(1:ix2,:,iz1:end,:);
    elseif strcmp(slice,'Left Leg')
        if strcmp(s(1),'b')
            z1 = 0.15;
            x1 = 0;
        elseif strcmp(s(1),'e')
            z1  = 0.1;
            x1 = 0;
        elseif strcmp(s(1),'d')
            z1  = 0.15;
            x1 = 0;
        end
        [~,iz1] = min(abs(r(1,1,:,3)-z1));
        [~,ix1] = min(abs(r(:,1,1,1)-x1));
        r= r(ix1:end,:,iz1:end,:);
        sigma_e = sigma_e(ix1:end,:,iz1:end,:);
        epsilon_r = epsilon_r(ix1:end,:,iz1:end,:);
    end
    idxS = find(abs(sigma_e(:)));
    xd = r(:,:,:,1);
    yd = r(:,:,:,2);
    x1  = min(xd(idxS));
    x2 = max(xd(idxS));
    y1  = min(yd(idxS));
    y2 = max(yd(idxS));
    [~,ix1] = min(abs(r(:,1,1,1)-x1));
    [~,ix2] = min(abs(r(:,1,1,1)-x2));
    [~,iy1] = min(abs(r(1,:,1,2)-y1));
    [~,iy2] = min(abs(r(1,:,1,2)-y2));
    r = r(ix1:ix2,iy1:iy2,:,:);
    epsilon_r = epsilon_r(ix1:ix2,iy1:iy2,:);
    sigma_e = sigma_e(ix1:ix2,iy1:iy2,:);
    
    [n1,n2,n3,~] = size(r);
    cent = r(round(n1/2),round(n2/2),round(n3/2),:);
    r(:,:,:,1) = r(:,:,:,1)-cent(1);
    r(:,:,:,2) = r(:,:,:,2)-cent(2);
    r(:,:,:,3) = r(:,:,:,3)-cent(3);
    
    if strcmp(s(1),'d')
        r(:,:,:,2)  = r(:,:,:,2) - 0.01;
        r(:,:,:,1)  = r(:,:,:,1) - 0.005;
    end
    rhbm_out.r = r;
    rhbm_out.epsilon_r = epsilon_r(end:-1:1,end:-1:1,end:-1:1);
    rhbm_out.sigma_e = sigma_e(end:-1:1,end:-1:1,end:-1:1);
    rho = zeros(size(rhbm_out.epsilon_r));
    rhbm_out.rho = rho;
    idxS = find(abs(rhbm_out.sigma_e(:)));
    rhbm_out.idxS = idxS;
    
end



