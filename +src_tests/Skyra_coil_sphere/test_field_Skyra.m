% MARIE 2.0 example

% addpath(genpath('../'));

% task_filename
task_path  = './+src_tests/Skyra_coil_sphere';
task_fname = 'Skyra_sphere.ini'; 


HFSS_files = {'SKYRA_shield_123MHz_Efield_TM_V4.fld',...
              'SKYRA_shield_123MHz_Hfield_TM_V4.fld'};


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

slice = 31;
freq  = 123e6; 

% %%
% 
% h_etot_axial    = simulator.show_total_electric_field('Axial', slice, 1, freq);
% h_etot_coronal  = simulator.show_total_electric_field('Coronal', slice, 1, freq);
% h_etot_sagittal = simulator.show_total_electric_field('Sagittal', slice, 1, freq);
% 
% %% show magnetic fields
% h_htot_axial    = simulator.show_total_magnetic_field('Axial', slice, 1, freq);
% h_htot_coronal  = simulator.show_total_magnetic_field('Coronal', slice, 1, freq);
% h_htot_sagittal = simulator.show_total_magnetic_field('Sagittal', slice, 1, freq);


axial     = @(Einc,slice)  sqrt(Einc(:,:,slice,1).^2 + Einc(:,:,slice,2).^2 + Einc(:,:,slice,3).^2);
coronal     = @(Einc,slice)  squeeze(sqrt(Einc(:,slice,:,1).^2 + Einc(:,slice,:,2).^2 + Einc(:,slice,:,3).^2));
sagittal     = @(Einc,slice)  squeeze(sqrt(Einc(slice,:,:,1).^2 + Einc(slice,:,:,2).^2 + Einc(slice,:,:,3).^2));

axial_lin    = @(E,cut) sqrt(squeeze(E(:,:,cut,1).^2) + squeeze(E(:,:,cut,5).^2) + squeeze(E(:,:,cut,9).^2));
coronal_lin  = @(E,cut) sqrt(squeeze(E(:,cut,:,1).^2) + squeeze(E(:,cut,:,5).^2) + squeeze(E(:,cut,:,9).^2));
sagittal_lin = @(E,cut) sqrt(squeeze(E(cut,:,:,1).^2) + squeeze(E(cut,:,:,5).^2) + squeeze(E(cut,:,:,9).^2));

%%
cut = slice;

H_tot = SOL_marie2.H_tot_;
E_tot = SOL_marie2.E_tot_;

H_marie = H_tot(:,:,:,:,1) + 1i * H_tot(:,:,:,:,2);
E_marie = E_tot(:,:,:,:,1) + 1i * E_tot(:,:,:,:,2);


PWX_flag = simulator.task_specs_.vie.PWX;

% mask out air voxels
if PWX_flag
    % mask out air voxels
    E_marie_axial = axial_lin(E_marie, cut);
    H_marie_axial = axial_lin(H_marie, cut);
    E_marie_coronal = coronal_lin(E_marie, cut);
    E_marie_sagittal = sagittal_lin(E_marie, cut);
else
    E_marie_cut = axial(E_marie, cut);
    E_marie_coronal = coronal(E_marie, cut);
    E_marie_sagittal = sagittal(E_marie, cut);
    H_marie_cut = axial(H_marie, cut);
end

idx_axial    = find(E_marie_axial.' == 0);
idx_coronal  = find(E_marie_coronal == 0); 
idx_sagittal = find(E_marie_sagittal == 0);

%% get electric and magnetic fields (HFSS)

[E_HFSS,~,~,~] = src_utils.read_fld_field(HFSS_files{1});
[H_HFSS,~,~,~] = src_utils.read_fld_field(HFSS_files{2});

% mask out air voxels
E_HFSS_axial = axial(E_HFSS, cut);
H_HFSS_axial = axial(H_HFSS, cut);
E_HFSS_coronal = coronal(E_HFSS, cut);
H_HFSS_coronal = coronal(H_HFSS, cut);
E_HFSS_sagittal = sagittal(E_HFSS, cut);
H_HFSS_sagittal = sagittal(H_HFSS, cut);

E_HFSS_axial(idx_axial) = 0;
H_HFSS_axial(idx_axial) = 0;
E_HFSS_coronal(idx_coronal) = 0;
H_HFSS_coronal(idx_coronal) = 0;
E_HFSS_sagittal(idx_sagittal) = 0;
H_HFSS_sagittal(idx_sagittal) = 0;
%% 

figure(); imagesc(abs(E_HFSS_sagittal));
caxis([0 1.4]);
colorbar;


figure(); imagesc(log10(abs(H_HFSS_coronal)));
colorbar;



