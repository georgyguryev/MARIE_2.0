function visualize_Network_param_matrix(E, freq_num, NP_title)

    % Visualizes the matrix of the S paramateres for a frequency point
    % Author: Ilias I. Giannakopoulos, Cambridge, MA, 2019
    
    
    % get scattering matrix for given frequency
    SPm = squeeze(E(:,:,freq_num));
        
%     MS = min(min(20*log10(abs(SPm))));
    
    % create a new figure
    figure
    
    % plot scattering matrix
    imagesc(20*log10(abs(SPm)));
%     for i = 1:size(SPm,1)
%         for j = 1:size(SPm,2)
%             text(i-0.5,j,num2str(20*log10(abs(SPm(i,j))), '%2.2f'), 'Color', [20*log10(abs(SPm(i,j)))/MS 1 1]);
%         end
%     end
    
    colorMap = [linspace(0,1,256)', zeros(256,2)];
    colormap(colorMap);
    colorbar;
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    title(NP_title)
    xlabel('Port #');
    ylabel('Port #');
    axis equal;
    
end