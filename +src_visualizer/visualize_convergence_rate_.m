function visualize_convergence_rate_(Cr,port_cut,freq)
    % Visualizes the convergence rate of the iterative solver
    % Author: Ilias I. Giannakopoulos, Cambridge, MA, 2019
    % Modified: Georgy Guryev, Cambridge, MA, 2019
    
    figure
    
    % get vector of residuals from the cell array
    res_vec = Cr{port_cut}; 
        
    semilogy(res_vec,'b');
    xlabel('Iteration Count');
    ylabel('Relative Residual');
    xlim([1, size(res_vec,1)]);
    grid on;
    axis tight;
    
end