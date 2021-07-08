function visualize_ConvergenceRate_(Cr,port_cut,freq_cut)

    % Visualizes the convergence rate of the iterative solver
    % Author: Ilias I. Giannakopoulos, Cambridge, MA, 2019
    
    figure
    
    if ndims(Cr) == 3
        Cr =  squeeze(Cr(:,port_cut,:));
        Cr =  squeeze(Cr(:,freq_cut));
    elseif ismatrix(Cr)
        Cr =  squeeze(Cr(:,port_cut));
    end
        
    semilogy(Cr,'b');
    xlabel('Iteration Count');
    ylabel('Relative Residue');
    grid on;
    axis tight;
    
end