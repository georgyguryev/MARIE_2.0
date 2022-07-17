% script is written to analyse and visualise results of simulations 

path_2_res  = './+src_tests/pFFT_vs_Basis/PWC/results/1e-3/';
path_2_figs = './+src_tests/pFFT_vs_Basis/PWC/results/figs/1e-3/';
data_files = {'Skyra_coil_RHBM_5mm_'};
          
figure_titles = {'Skyra coil with RHBM HT 5mm '};
          
iter_methods = {'GMRES_DR.mat',  'TFQMR.mat'};

figure_it_methods = {' GMRES DR',' TFQMR'};

          
%% read files

% loop over data root names
for i = 1:length(data_files)
    
    res_root_name = data_files{i};
    
    
    % loop over iterative methods
    for j = 1:length(iter_methods)
        
        iter_method   = iter_methods{j}; 
        
        % get full file name 
        res_file_name = strcat(res_root_name, iter_method);
        fig_file_name = strcat(figure_titles{i}, figure_it_methods{j});
        time_fig_name = strcat(strcat('Runtime ', figure_titles{i}), figure_it_methods{j});
        
        % if file with this name exists run analysis
        if 2 == exist(res_file_name, 'file')
            
            res_fig = figure();
            
            % load data
            load(res_file_name);
            
            Sol_exp_coup   = SOL_marie2{1};
            Sol_basis_coup = SOL_marie2{2}; 
            Sol_exp_decoup = SOL_marie2{3};
            Sol_basis_decoup = SOL_marie2{4};
            
            % get residual vectors
            res_vec_exp_coup  = Sol_exp_coup.res_vec_;
            res_vec_basis_coup = Sol_basis_coup.res_vec_;
            res_vec_exp_decoup = Sol_exp_decoup.res_vec_;
            res_vec_basis_decoup = Sol_exp_decoup.res_vec_;
            
            % get polarization currents
            Jb_exp_coup   = Sol_exp_coup.Jb_plrz_vec_;
            Jb_basis_coup = Sol_basis_coup.Jb_plrz_vec_;
            Jb_exp_decoup = Sol_exp_decoup.Jb_plrz_vec_;
            Jb_basis_decoup = Sol_basis_decoup.Jb_plrz_vec_;

            
            semilogy(res_vec_exp_coup{1}, '-x');
            hold on;
            semilogy(res_vec_basis_coup{1}, '-.');
            semilogy(res_vec_exp_decoup{1}, '-x');
            semilogy(res_vec_basis_decoup{1}, '-.');
            
            legend_exp_coup = sprintf('Explicit Dense method');
            legend_basis_coup = sprintf('Explicit Basis method');
            legend_exp_decoup = sprintf('Implicit Dense method');
            legend_basis_decoup = sprintf('Implicit Basis method');

            legend(legend_exp_coup, legend_basis_coup,...
                   legend_exp_decoup, legend_basis_decoup);
               
            title(fig_file_name);
            xlabel('Number of iterations');
            ylabel('relative residual');
            hold off;
            
            fig_full_fname = fullfile(path_2_figs, strcat(fig_file_name, '.fig'));
            png_full_fname = fullfile(path_2_figs, strcat(fig_file_name, '.png'));
            
            % save figures as .fig and .png
            savefig(res_fig, fig_full_fname);
            saveas(res_fig,png_full_fname);
            
            close all;
            
            % visualize average runtime
            t_average     = [Sol_exp_coup.average_runtime();...
                Sol_basis_coup.average_runtime();...
                Sol_exp_decoup.average_runtime();...
                Sol_basis_decoup.average_runtime()];
            
            time_fig = figure(); 
            
            bar([1], t_average);
            set(gca, 'yscale', 'log');  
            hold on;
            xlim = get(gca, 'xlim');
            ylim([0.1, 2*round(max(t_average))]);
            plot(xlim, [min(t_average) min(t_average)], '-.k');
            ylabel('Average runtime, sec');
            xlabel('Solver configurations');
            legend(legend_exp_coup,legend_basis_coup,...
                legend_exp_decoup, legend_basis_decoup, 'Best runtime');
            
            title(time_fig_name);
            
            time_fig_full_fname = fullfile(path_2_figs, strcat(time_fig_name, '.fig'));
            time_png_full_fname = fullfile(path_2_figs, strcat(time_fig_name, '.png'));
            
            % save figures as .fig and .png
            savefig(time_fig, time_fig_full_fname);
            saveas(time_fig, time_png_full_fname);
            
            close all;
        end
        
        % visualize average iteration count
        
        
        
    end
    
end