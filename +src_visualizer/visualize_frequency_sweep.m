function visualize_frequency_sweep(Sparam,i,j,freqs)

    % Visualizes one S paramater S(i,j) for all the frequency points
    % Author: Ilias I. Giannakopoulos, Cambridge, MA, 2019

    figure
    
    % select entry of interest
    SPm = squeeze(Sparam(i,j,:));
    
    % convert S parameters to dBs
    Splot = 20*log10(abs(SPm));
    
    % plot the frequency sweep
    plot(freqs,Splot,'b');
    xlabel('frequencies');
    ylabel('S (dB)');
    grid on;
end