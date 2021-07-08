function visualize_properties_(EP,sac,cut)

    % Visualizes a slice of one electrical property of the body model
    % Author: Ilias I. Giannakopoulos, Cambridge, MA, 2019
    
    % check if data is empty
    if isempty(EP)
        warning("Can not visualize electrical properties: Tensor is empty! ");
        return
    end
    
    figure();
    
    % crop slice
    if strcmp(sac,'Sagittal')
        
        EP = squeeze(EP(cut,:,:));
        
    elseif strcmp(sac,'Coronal')
        
        EP = squeeze(EP(:,cut,:));
        
    elseif strcmp(sac,'Axial')
        
        EP = squeeze(EP(:,:,cut));
        
    end
    
    imagesc(abs(EP));
    colormap hot;
    colorbar
    axis equal;
end