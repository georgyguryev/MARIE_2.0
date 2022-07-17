% this script is used to compare frequency sweep of S parameters as a
% function of Birdcage coil refinement 

path_to_files = fullfile('.', '+src_tests', 'pFFT_vs_Basis', 'PWC', 'sweep', '1e-3');
 
% file_coarse = 'Bcage_Cylinder_0.145R_0.3H_res_5mm_GMRES_DR_coarse.mat';
% file_normal = 'Bcage_Cylinder_0.145R_0.3H_res_5mm_GMRES_DR.mat';
% file_refined = 'Bcage_Cylinder_0.145R_0.3H_res_5mm_TFQMR.mat';


file_coarse = 'Bcage_RHBM_res_5mm_GMRES_DR_coarse.mat';
file_normal = 'Bcage_RHBM_res_5mm_GMRES_DR.mat';
file_refined = 'Bcage_RHBM_res_5mm_TFQMR.mat';


%%

fname_coarse = fullfile(path_to_files, file_coarse);
fname_normal = fullfile(path_to_files, file_normal);
fname_refined = fullfile(path_to_files, file_refined);


%% load files 

load(fname_coarse);
Sol_coarse = SOL_marie2{1};
clear SOL_marie2;

load(fname_normal);
Sol_normal = SOL_marie2{1};
clear SOL_marie2;


load(fname_refined);
Sol_refined = SOL_marie2{1};
clear SOL_marie2;


%% 

freq = 1e6 * linspace(275, 325, 21);

Sij_coil = @(Sol, i,j) squeeze(Sol.Network_param_.coil.S(i,j,:));
Sij_src = @(Sol, i,j) squeeze(Sol.Network_param_.src.S(i,j,:));

legend_str = {'coarse mesh', 'normal mesh', 'refined mesh'};
%%


for i=1:2
    for j = 1:2

        title_str = sprintf('S_{%d %d}, dB', i,j);
        
        figure(); plot(freq, 20 * log10(abs(Sij_coil(Sol_coarse,i,j))));
        hold on;
        plot(freq, 20 * log10(abs(Sij_coil(Sol_normal,i,j))), '-o');
        plot(freq, 20 * log10(abs(Sij_coil(Sol_refined,i,j))), '-x');
        
        legend(legend_str{:})
        xlabel('Frequency, Hz');
        ylabel('S param, dB');
        grid on;
        title(title_str);

    end
end

