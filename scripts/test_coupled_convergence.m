% This script is designed to compare the convergence properties of the
% coupled and decoupled basis-based solvers 

% provide full name and path to the ini file 
ini_filename = './scripts/test_coupled_convergence_coupled.ini'; 

% create simulator
simulator = src_task.EM_simulator();

% load simulation specs 
simulator.load_task(ini_filename);
%% visualize problem to be solved

% simulator.show_coil_and_body();
% 
% keyboard; 

%% run simulation

simulator.run_task();

%% get solution

Sol = simulator.get_solution();

%%
slice = 52;
% freq  = simulator.; 
res = 1;

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

% h_etot_axial    = simulator.show_total_electric_field('Axial', slice, 1, freq);
% title(strcat(title_Etot, title_axial), 'Interpreter', 'latex');
h_etot_coronal  = simulator.show_total_electric_field('Coronal', slice, 1, freq);
title(strcat(title_Etot, title_coronal), 'Interpreter', 'latex');
% h_etot_sagittal = simulator.show_total_electric_field('Sagittal', slice,
% 1, freq); title(strcat(title_Etot, title_sagittal), 'Interpreter',
% 'latex');

%% show magnetic fields
% h_htot_axial    = simulator.show_total_magnetic_field('Axial', slice, 1, freq);
% title(strcat(title_Htot, title_axial), 'Interpreter', 'latex');
h_htot_coronal  = simulator.show_total_magnetic_field('Coronal', slice, 1, freq);
title(strcat(title_Htot, title_coronal), 'Interpreter', 'latex');
% h_htot_sagittal = simulator.show_total_magnetic_field('Sagittal', slice, 1, freq);
% title(strcat(title_Htot, title_sagittal), 'Interpreter', 'latex');


