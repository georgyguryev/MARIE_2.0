function h = visualize_fields_(E,sac,cut,port_cut,freq_cut)

    % Visualizes the slice of an electromagnetic field measurement magnitude
    % Author: Ilias I. Giannakopoulos, Cambridge, MA, 2019
    
    h = figure();
    
    if ndims(E) > 4
        % crop port and frequency points
        E =  squeeze(E(:,:,:,:,port_cut,freq_cut));
        
    end
    
    % calculate magnitude if needed for PWC, or PWL measurements
    if ndims(E) > 3
        
        [~,~,~,n] = size(E);
        
        if n == 3
        
            E = sqrt(abs(E(:,:,:,1)).^2+abs(E(:,:,:,2)).^2+abs(E(:,:,:,3)).^2);

        elseif n == 12
        
            E = sqrt(abs(E(:,:,:,1)).^2+abs(E(:,:,:,5)).^2+abs(E(:,:,:,9)).^2);
        
        end
        
    end
    
    sac = lower(sac);
    
    % crop slice
    if strcmp(sac,'sagittal')
        
        H = squeeze(E(cut,:,:));
        
    elseif strcmp(sac,'coronal')
        
        H = squeeze(E(:,cut,:));
        
    elseif strcmp(sac,'axial')
        
        H = squeeze(E(:,:,cut));
        
    end
    
    imagesc(H);
    colormap hot;
    colorbar
    axis equal;
    
end