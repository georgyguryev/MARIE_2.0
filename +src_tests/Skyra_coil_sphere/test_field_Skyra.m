% MARIE 2.0 example

% addpath(genpath('../'));

% task_filename
task_path  = './+src_tests/Skyra_coil_sphere';
task_fname = 'Skyra_sphere.ini'; 

% get full file name with a path
task_specs = fullfile(task_path, task_fname);

%% Load task

% create simulator
simulator = src_task.EM_simulator();

% construct appropriate task from task_filename with a factory method
simulator.load_task(task_specs);

%% Run task

% profile on;

% run simulation task
simulator.run_task();

%% Return results of simulation

SOL_marie2 = simulator.get_solution();
 
%% Run tests for implants 

slice = 31;
freq  = 123e6; 

%%

h_etot_axial    = simulator.show_total_electric_field('Axial', slice, 1, freq);
h_etot_coronal  = simulator.show_total_electric_field('Coronal', slice, 1, freq);
h_etot_sagittal = simulator.show_total_electric_field('Sagittal', slice, 1, freq);

%% show magnetic fields
h_htot_axial    = simulator.show_total_magnetic_field('Axial', slice, 1, freq);
h_htot_coronal  = simulator.show_total_magnetic_field('Coronal', slice, 1, freq);
h_htot_sagittal = simulator.show_total_magnetic_field('Sagittal', slice, 1, freq);


