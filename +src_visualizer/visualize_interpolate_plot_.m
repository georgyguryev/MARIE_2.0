function fig_handle = visualize_interpolate_plot_(E,sac,cut,port_cut,freq_cut,stepi)

    % Performs an interpolation plot for a slice of an electromagnetic field measurement magnitude
    % Author: Ilias I. Giannakopoulos, Cambridge, MA, 2019
    
    fig_handle = figure();
    
    % crop port and frequency points
    E =  squeeze(E(:,:,:,:,port_cut,freq_cut));
    
    [L,M,N,~] = size(E);
    
    L_f = zeros(L,stepi);
    M_f = zeros(M,stepi);
    N_f = zeros(N,stepi);
    
    for i = 1:stepi
        L_f(:,i) = i:stepi:stepi*L-(stepi-i);
        M_f(:,i) = i:stepi:stepi*M-(stepi-i);
        N_f(:,i) = i:stepi:stepi*N-(stepi-i);
    end
    
    E_f = zeros(stepi*L,stepi*M,stepi*N,3);
    
    for l = 1:stepi
        for m = 1:stepi
            for n = 1:stepi
                
                l_step = -(stepi-(2*l-1))/(2*stepi);
                m_step = -(stepi-(2*m-1))/(2*stepi);
                n_step = -(stepi-(2*n-1))/(2*stepi);
                
                for i = 1:3   
                    
                    E_f (L_f(:,l),M_f(:,m),N_f(:,n),i) = E(:,:,:,(i-1)*4+1) / stepi^3 + ...
                                                          l_step * E(:,:,:,(i-1)*4+2) + ...
                                                          m_step * E(:,:,:,(i-1)*4+3) + ...
                                                          n_step *E(:,:,:,(i-1)*4+4);
                    
                end
                
            end
        end
    end
    
    % calculate magnitude of PWL measurements
    E2 = sqrt(abs(E_f(:,:,:,1)).^2+abs(E_f(:,:,:,2)).^2+abs(E_f(:,:,:,3)).^2);
    
    % crop slice
    if strcmp(sac,'Sagittal')
    
        cut = cut*stepi;
        if cut>L
            cut = L;
        end
        
        H = squeeze(E2(cut,:,:));
        
    elseif strcmp(sac,'Coronal')
        
        cut = cut*stepi;
        if cut > M
            cut = N;
        end
        
        H = squeeze(E2(:,cut,:));
        
    elseif strcmp(sac,'Axial')
        
        cut = cut*stepi;
        if cut>N
            cut = N;
        end
        
        H = squeeze(E2(:,:,cut));
        
    end
    
    imagesc(H);
    colormap hot;
    colorbar
    axis equal;
end