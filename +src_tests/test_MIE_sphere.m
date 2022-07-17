% close all;


res = 0.0025;
slice = 40; 

% set paths and filenames for MIE and MARIE 2.0
path_to_file = fullfile('data', 'bodies');
body_fname   = sprintf('Sphere_rad_0.100_res_%.4f.mat', res);
body_full_filename = fullfile(path_to_file, body_fname);

%% prepare MARIE 2.0

% task_filename
task_path  = './+src_tests';
figs_path  = './+src_tests/Mie_validation_figs/';
png_path   = fullfile(figs_path, 'png');
fig_path   = fullfile(figs_path, 'fig');
task_fname = 'VIE_MIE_validation.ini';

% get full file name with a path
task_specs = fullfile(task_path, task_fname);

%%
freq = 298e6;
sphere_rad = 0.10-eps;

src_utils.EM_utils(freq);


%% load file 
load(body_full_filename);

[L,M,N,~] = size(r);

i_c = sub2ind([L,M,N], round(L/2), round(M/2), round(N/2));

er_c = epsilon_r(i_c);
se_c = sigma_e(i_c);


%% run MIE series code 

profile on;

[E_Mie, H_Mie] = MIE_SERIES(r, sphere_rad, er_c, se_c,freq);

profile off;
profile viewer;

%%
axial    = @(E,cut) sqrt(squeeze(E(:,:,cut,1).^2) + squeeze(E(:,:,cut,2).^2) + squeeze(E(:,:,cut,3).^2));
coronal  = @(E,cut) sqrt(squeeze(E(:,cut,:,1).^2) + squeeze(E(:,cut,:,2).^2) + squeeze(E(:,cut,:,3).^2));
sagittal = @(E,cut) sqrt(squeeze(E(cut,:,:,1).^2) + squeeze(E(cut,:,:,2).^2) + squeeze(E(cut,:,:,3).^2));


axial_lin    = @(E,cut) sqrt(squeeze(E(:,:,cut,1).^2) + squeeze(E(:,:,cut,5).^2) + squeeze(E(:,:,cut,9).^2));
coronal_lin  = @(E,cut) sqrt(squeeze(E(:,cut,:,1).^2) + squeeze(E(:,cut,:,5).^2) + squeeze(E(:,cut,:,9).^2));
sagittal_lin = @(E,cut) sqrt(squeeze(E(cut,:,:,1).^2) + squeeze(E(cut,:,:,5).^2) + squeeze(E(cut,:,:,9).^2));


%%

e_mie_ax  = figure(); imagesc(20 * log10(abs(axial(E_Mie, slice))));
e_mie_cor = figure(); imagesc(20 * log10(abs(coronal(E_Mie, slice))));
e_mie_sag = figure(); imagesc(20 * log10(abs(sagittal(E_Mie, slice))));

h_mie_ax  = figure();  imagesc(20 * log10(abs(axial(H_Mie, slice))));
h_mie_cor = figure(); imagesc(20 * log10(abs(coronal(H_Mie, slice))));
h_mie_sag = figure(); imagesc(20 * log10(abs(sagittal(H_Mie, slice))));



%% 

str_e_axial = sprintf('Etot_Mie_Sphere_0.1R_%.4fres_axial_%d.', res,slice);
str_e_coron = sprintf('Etot_Mie_Sphere_0.1R_%.4fres_coron_%d.', res, slice);
str_e_sagit = sprintf('Etot_Mie_Sphere_0.1R_%.4fres_sagit_%d.', res, slice);

str_h_axial = sprintf('Htot_Mie_Sphere_0.1R_%.4fres_axial_%d.', res,slice);
str_h_coron = sprintf('Htot_Mie_Sphere_0.1R_%.4fres_coron_%d.', res, slice);
str_h_sagit = sprintf('Htot_Mie_Sphere_0.1R_%.4fres_sagit_%d.', res, slice);


saveas(e_mie_ax,  fullfile(png_path, strcat(str_e_axial, 'png')));
saveas(e_mie_cor, fullfile(png_path, strcat(str_e_coron, 'png')));
saveas(e_mie_sag, fullfile(png_path, strcat(str_e_sagit, 'png')));

saveas(h_mie_ax,  fullfile(png_path, strcat(str_h_axial, 'png')));
saveas(h_mie_cor, fullfile(png_path, strcat(str_h_coron, 'png')));
saveas(h_mie_sag, fullfile(png_path, strcat(str_h_sagit, 'png')));


%% save as fig

savefig(e_mie_ax,  fullfile(fig_path, strcat(str_e_axial, 'fig')));
savefig(e_mie_cor, fullfile(fig_path, strcat(str_e_coron, 'fig')));
savefig(e_mie_sag, fullfile(fig_path, strcat(str_e_sagit, 'fig')));

savefig(h_mie_ax,  fullfile(fig_path, strcat(str_h_axial, 'fig')));
savefig(h_mie_cor, fullfile(fig_path, strcat(str_h_coron, 'fig')));
savefig(h_mie_sag, fullfile(fig_path, strcat(str_h_sagit, 'fig')));


close all;

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

E_tot = SOL_marie2.E_tot_;
H_tot = SOL_marie2.H_tot_;


