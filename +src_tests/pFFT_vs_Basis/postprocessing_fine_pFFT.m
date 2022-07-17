% script is written to analyse and visualise results of simulations 

path_pwl_res  = './+src_tests/pFFT_vs_Basis/PWL/results/1e-3/';
path_pwc_res  = './+src_tests/pFFT_vs_Basis/PWC/results/1e-3/'; 

path_2_figs = './+src_tests/pFFT_vs_Basis/PWL/results/figs/1e-3/';
data_files = {'Bcage_Cylinder_0.1R_0.3H_res_5mm_',... 
              'Bcage_Cylinder_0.145R_0.3H_res_5mm_',...
              'Bcage_RHBM_res_5mm_',...
              'Triangular_coil_Cylinder_0.1R_0.3H_res_5mm_',...
              'Triangular_coil_RHBM_res_5mm_'};
          
figure_titles = {'pFFT Fine Bcage with Cylinder (rad = 0.1, res = 5mm)  ',...
                 'pFFT Fine Bcage with Cylinder (rad = 0.145, res = 5mm)  ',...
                 'pFFT Fine Bcage with RHBM (res = 5mm)  ',...
                 'pFFT Triangular coil with Cylinder (rad = 0.1, res = 0.005)',...
                 'pFFT Triangular coil with RHBM (res = 0.005)'};
          
iter_methods = {'GMRES_DR.mat',  'TFQMR.mat'};
fine_iter_methods = {'GMRES_DR_pFFT.mat',  'TFQMR_pFFT.mat'};

figure_it_methods = {'GMRES DR','TFQMR'};

%%

% loop over data root names
for i = 1:length(data_files)
    
    res_root_name = data_files{i};
    
    
    % loop over iterative methods
    for j = 1:length(iter_methods)
        
        iter_method_pwl = iter_methods{j}; 
        iter_method_pwc = fine_iter_methods{j};
        
        % get full file name 
        pwl_file_name = fullfile(path_pwl_res,  strcat(res_root_name, iter_method_pwl));
        pwc_file_name = fullfile(path_pwc_res, strcat(res_root_name, iter_method_pwc));
        fig_file_name = strcat(figure_titles{i}, figure_it_methods{j});
        time_fig_name = strcat(strcat('Runtime ', figure_titles{i}), figure_it_methods{j});
        
        % if file with this name exists run analysis
        if 2 == exist(pwl_file_name, 'file')
            
            res_fig = figure();
            
            % load data
            load(pwl_file_name);
            
            Sol_pFFT_exp_pwl  = SOL_marie2{1};
            Sol_pFFT_imp_pwl = SOL_marie2{2};
            
            % get residual vectors
            
            res_vec_pFFT_exp_pwl = Sol_pFFT_exp_pwl.res_vec_;
            res_vec_pFFT_imp_pwl = Sol_pFFT_imp_pwl.res_vec_;
            
            % get polarization currents
            Jb_pFFT_exp_pwl  = Sol_pFFT_exp_pwl.Jb_plrz_vec_;
            Jb_pFFT_imp_pwl  = Sol_pFFT_imp_pwl.Jb_plrz_vec_;
            
            load(pwc_file_name);
            
            Sol_pFFT_exp_pwc = SOL_marie2{1};
            Sol_pFFT_imp_pwc = SOL_marie2{2};
            
            res_vec_pFFT_exp_pwc = Sol_pFFT_exp_pwc.res_vec_;
            res_vec_pFFT_imp_pwc = Sol_pFFT_imp_pwc.res_vec_;
            
            % get polarization currents
            Jb_pFFT_exp = Sol_pFFT_exp_pwc.Jb_plrz_vec_;
            Jb_pFFT_imp = Sol_pFFT_imp_pwc.Jb_plrz_vec_;
                                 

            semilogy(res_vec_pFFT_exp_pwl{1}, '-o');
            hold on;
            semilogy(res_vec_pFFT_imp_pwl{1}, '-o');
            semilogy(res_vec_pFFT_exp_pwc{1}, '-x');
            semilogy(res_vec_pFFT_imp_pwc{1}, '-x');
            

            legend_pFFT_exp_pwl = sprintf('Explicit pFFT method, pwl');
            legend_pFFT_imp_pwl = sprintf('Implicit pFFT method, pwl');
            legend_pFFT_exp_pwc = sprintf('Explicit pFFT method, pwc');
            legend_pFFT_imp_pwc = sprintf('Implicit pFFT method, pwc');
            
            legend(legend_pFFT_exp_pwl, legend_pFFT_imp_pwl,...
                   legend_pFFT_exp_pwc, legend_pFFT_imp_pwc);   
               
            title(fig_file_name);
            hold off;
            
            fig_full_fname = fullfile(path_2_figs, strcat(fig_file_name, '.fig'));
            png_full_fname = fullfile(path_2_figs, strcat(fig_file_name, '.png'));
            
            % save figures as .fig and .png
            savefig(res_fig, fig_full_fname);
            saveas(res_fig,png_full_fname);
            
            close all;
            
            % visualize average runtime
            t_average     = [Sol_pFFT_exp_pwl.average_runtime();...
                             Sol_pFFT_imp_pwl.average_runtime();...
                             Sol_pFFT_exp_pwc.average_runtime();...
                             Sol_pFFT_imp_pwc.average_runtime()];
            
            time_fig = figure(); 
            bar([1], t_average);
%             set(gca, 'yscale', 'log');  

            hold on;
            xlim = get(gca, 'xlim');
            plot(xlim, [min(t_average) min(t_average)], '-.k');
            ylim([1, 2*round(max(t_average))]);
            ylabel('Average runtime, sec');
            xlabel('Solver configurations');
            
            legend(legend_pFFT_exp_pwl,legend_pFFT_imp_pwl,...
                   legend_pFFT_exp_pwc, legend_pFFT_imp_pwc,...
                   'Best runtime');
            
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
