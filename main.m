% MARIE 2.0 example

% addpath(genpath('../'));

% task_filename
task_path  = './+src_tests';
task_fname = 'sphere_small_pwc_lr.ini'; 

path_to_figs = '/home/georgy/ISMRM_2022/figures/fig/Billie_pTx_8ch';

% get full file name with a path
task_specs = fullfile(task_path, task_fname);

%% Load task

% create simulator
simulator = src_task.EM_simulator();

% construct appropriate task from task_filename with a factory method
simulator.load_task(task_specs);

%% Run task

% run simulation task
simulator.run_task();

%% Return results of simulation

SOL_marie2 = simulator.get_solution();
 
%% Run tests for implants 

slice = 36;
freq  = 400.2e6; 
res = 5;

%%
title_Etot  = '$${\bf E_{tot}}, [V/m]$$';
title_Htot  = '$${\bf H_{tot}}, [A/m]$$';
title_axial    = sprintf('\nAxial cut %d, %dmm res', slice, res); 
title_coronal  = sprintf('\nCoronal cut %d, %dmm res', slice, res);
title_sagittal = sprintf('\nSagittal cut %d, %dmm res', slice, res); 


figfile_Etot = 'E_tot_';
figfile_Htot = 'H_tot_';
figfile_axial    = sprintf('axial_cut_%d_res_%dmm', slice, res);
figfile_coronal  = sprintf('coronal_cut_%d_res_%dmm', slice, res);
figfile_sagittal = sprintf('sagittal_cut_%d_res_%dmm', slice, res);


%% show electric fields 

h_etot_axial    = simulator.show_total_electric_field('Axial', slice, 1, freq);
title(strcat(title_Etot, title_axial), 'Interpreter', 'latex');
h_etot_coronal  = simulator.show_total_electric_field('Coronal', slice, 1, freq);
title(strcat(title_Etot, title_coronal), 'Interpreter', 'latex');
h_etot_sagittal = simulator.show_total_electric_field('Sagittal', slice, 1, freq);
title(strcat(title_Etot, title_sagittal), 'Interpreter', 'latex');

%% show magnetic fields

h_htot_axial    = simulator.show_total_magnetic_field('Axial', slice, 1, freq);
title(strcat(title_Htot, title_axial), 'Interpreter', 'latex');
h_htot_coronal  = simulator.show_total_magnetic_field('Coronal', slice, 1, freq);
title(strcat(title_Htot, title_coronal), 'Interpreter', 'latex');
h_htot_sagittal = simulator.show_total_magnetic_field('Sagittal', slice, 1, freq);
title(strcat(title_Htot, title_sagittal), 'Interpreter', 'latex');

% %%  electric fields 
% savefig(h_etot_axial, fullfile(path_to_figs, strcat(figfile_Etot, figfile_axial)));
% savefig(h_etot_coronal, fullfile(path_to_figs, strcat(figfile_Etot, figfile_coronal)));
% savefig(h_etot_sagittal, fullfile(path_to_figs, strcat(figfile_Etot, figfile_sagittal)));
% 
% %% magnetic fields
% savefig(h_htot_axial, fullfile(path_to_figs, strcat(figfile_Htot, figfile_axial)));
% savefig(h_htot_coronal, fullfile(path_to_figs, strcat(figfile_Htot, figfile_coronal)));
% savefig(h_htot_sagittal, fullfile(path_to_figs, strcat(figfile_Htot, figfile_sagittal)));



