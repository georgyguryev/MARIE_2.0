% MARIE 2.0 example

% addpath(genpath('../'));

% close all;

% task_filename
task_path  = './+src_tests/16ch_coil_loop1_tuned_matched';
task_fname = '16ch_loop1_decoupled_gaps_curved_V1.ini'; 

HFSS_files = {'Efield_HFSS_16ch_loop1_decoupled_gaps_50mmsphere_2p5mmres.fld',...
              'Hfield_HFSS_16ch_loop1_decoupled_gaps_50mmsphere_2p5mmres.fld'};

% get full file name with a path
task_specs = fullfile(task_path, task_fname);

%%

cut = 21;

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

%% set axial 

axial     = @(Einc,slice)  sqrt(Einc(:,:,slice,1).^2 + Einc(:,:,slice,2).^2 + Einc(:,:,slice,3).^2);
axial_lin = @(Einc, slice) sqrt(Einc(:,:,slice,1).^2 + Einc(:,:,slice,5).^2 + Einc(:,:,slice,9).^2);

%% get electric and magnetic fields (MARIE 2.0)

H_marie = SOL_marie2.H_tot_;
E_marie = SOL_marie2.E_tot_;


% mask out air voxels
E_marie_cut = axial(E_marie, cut);
H_marie_cut = axial(H_marie, cut);

idx = find(E_marie_cut == 0);

%% get electric and magnetic fields (HFSS)

[E_HFSS,~,~,~] = src_utils.read_fld_field(HFSS_files{1});
[H_HFSS,~,~,~] = src_utils.read_fld_field(HFSS_files{2});

% mask out air voxels
E_HFSS_cut = axial(E_HFSS, cut);
H_HFSS_cut = axial(H_HFSS, cut);

E_HFSS_cut(idx) = 0;
H_HFSS_cut(idx) = 0;

%% 

figure(); imagesc(abs(E_HFSS_cut));
caxis([0 2.1]);
colorbar;

figure(); imagesc(abs(E_marie_cut));
caxis([0 2.1]);
colorbar; 

figure(); imagesc(abs(H_HFSS_cut));
colorbar;

figure(); imagesc(abs(H_marie_cut));
colorbar;
