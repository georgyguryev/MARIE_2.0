function [rhbm_out] = RHBM_RHBM_pad_(s,slice,pos,er,se)
    
    % Crops a Virtual family model and generates a realistic human body
    % model by attaching a hmogeneous dielectric pad on one side of the
    % head
    % Author: Ilias I. Giannakopoulos, Cambridge, MA, 2019
    
    rhbm_out.name = strcat(s(1:end-4),'_',slice);
    load(strcat('./src/data/bodies/Virtual_Family/',s),'epsilon_r','sigma_e','r');
    
    %% Cut the air around the whole body
    yd = r(:,:,:,2);
    zd = r(:,:,:,3);
    idxS = find(abs(sigma_e(:)));
    y1  = min(yd(idxS));
    y2 = max(yd(idxS));
    z1  = min(zd(idxS));
    z2 = max(zd(idxS));
    [~,iy1] = min(abs(r(1,:,1,2)-y1));
    [~,iy2] = min(abs(r(1,:,1,2)-y2));
    [~,iz1] = min(abs(r(1,1,:,3)-z1));
    [~,iz2] = min(abs(r(1,1,:,3)-z2));
    
    %% Take the head
    r = r(:,iy1:iy2,iz1:iz2,:);
    epsilon_r = epsilon_r(:,iy1:iy2,iz1:iz2);
    sigma_e = sigma_e(:,iy1:iy2,iz1:iz2);
    idxS = find(abs(sigma_e(:)));

    if strcmp(s(1),'b')
        z2  = -0.52;
        y_cut_2 = 0.04;
        y_cut_1 = 0.04;
    elseif strcmp(s(1),'e')
        z2  = -0.58;
        y_cut_2 = 0.04;
        y_cut_1 = 0.04;
    elseif strcmp(s(1),'d')
        z2  = -0.67;
        y_cut_2 = 0.08;
        y_cut_1 = 0.04;
    end
    y1  = min(yd(idxS));
    y2 = max(yd(idxS));
    [~,iy1] = min(abs(r(1,:,1,2)-y1));
    [~,iy2] = min(abs(r(1,:,1,2)-y2));
    [~,iz2] = min(abs(r(1,1,:,3)-z2));
    r= r(:,iy1:iy2,1:iz2,:);
    sigma_e = sigma_e(:,iy1:iy2,1:iz2);
    epsilon_r = epsilon_r(:,iy1:iy2,1:iz2);
    idxS = find(abs(sigma_e(:)));
   
    [n1,n2,n3,~] = size(r);
    z_cut = 0.02;
    x_cut = 0.01;
    xd = r(:,:,:,1);
    yd = r(:,:,:,2);
    zd = r(:,:,:,3);
    y1  = min(yd(idxS));
    y2 = max(yd(idxS));
    z1  = min(zd(idxS));
    z2 = max(zd(idxS));
    mask = zeros(n1,n2,n3);
    if strcmp(pos,'Right')        
        for i = 1:n1
            for j = 1:n2
                for k = 1:n3
                    if (yd(i,j,k) < y2-y_cut_2) && (yd(i,j,k) > y1+y_cut_1) && (zd(i,j,k) < z2-z_cut) && (zd(i,j,k) > z1+z_cut)
                        xdi = xd(:,j,k);
                        sigma_ei = sigma_e(:,j,k);
                        idxSx = find(abs(sigma_ei(:)));
                        x_bound  = max(xdi(idxSx));
                        if isempty(x_bound)
                            continue;
                        end
                        if (xd(i,j,k) > x_bound) && (xd(i,j,k) < x_bound+x_cut)
                            mask(i,j,k) = 1;
                        end
                    end
                end
            end
        end
    elseif strcmp(pos,'Left')
        for i = 1:n1
            for j = 1:n2
                for k = 1:n3
                    if (yd(i,j,k) < y2-y_cut_2) && (yd(i,j,k) > y1+y_cut_1) && (zd(i,j,k) < z2-z_cut) && (zd(i,j,k) > z1+z_cut)
                        xdi = xd(:,j,k);
                        sigma_ei = sigma_e(:,j,k);
                        idxSx = find(abs(sigma_ei(:)));
                        x_bound  = min(xdi(idxSx));
                        if isempty(x_bound)
                            continue;
                        end
                        if (xd(i,j,k) < x_bound) && (xd(i,j,k) > x_bound-x_cut)
                            mask(i,j,k) = 1;
                        end
                    end
                end
            end
        end
    end
    sigma_e = sigma_e+mask*se;
    epsilon_r = epsilon_r+mask*er;
    
    %% Cut the extra air
    idxS = find(abs(sigma_e(:)));
    xd = r(:,:,:,1);
    x1  = min(xd(idxS));
    x2 = max(xd(idxS));
    [~,ix1] = min(abs(r(:,1,1,1)-x1));
    [~,ix2] = min(abs(r(:,1,1,1)-x2));
    r = r(ix1:ix2,:,:,:);
    epsilon_r = epsilon_r(ix1:ix2,:,:);
    sigma_e = sigma_e(ix1:ix2,:,:);
    mask = mask(ix1:ix2,:,:);
    
    %% Centering
    [n1,n2,n3,~] = size(r);
    cent = r(round(n1/2),round(n2/2),round(n3/2),:);
    r(:,:,:,1) = r(:,:,:,1)-cent(1);
    r(:,:,:,2) = r(:,:,:,2)-cent(2);
    r(:,:,:,3) = r(:,:,:,3)-cent(3);
    
    %% Output
    rhbm_out.r = r;
    rhbm_out.epsilon_r = epsilon_r(end:-1:1,end:-1:1,end:-1:1);
    rhbm_out.sigma_e = sigma_e(end:-1:1,end:-1:1,end:-1:1);
    rho = zeros(size(rhbm_out.epsilon_r));
    rhbm_out.rho = rho;
    idxS = find(abs(rhbm_out.sigma_e(:)));
    rhbm_out.idxS = idxS;
    mask = mask(end:-1:1,end:-1:1,end:-1:1);
    idxS2 = find(mask(:));
    rhbm_out.idxS2 = idxS2;
end



