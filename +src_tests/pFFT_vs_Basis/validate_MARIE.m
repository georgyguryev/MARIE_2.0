close all;

path_to_hffs    = fullfile('~','Documents','MIT', 'PhD', 'projects', 'MARIE_v2.0', 'Eugene_data',...
                           'MARIE_comparisons','MARIE_comparisons', 'Sphere_15cm_rad_5mm_res');
path_to_results = fullfile('+src_tests', 'pFFT_vs_Basis', 'PWC', 'results', '1e-3');

file_name       = 'Skyra_coil_RHBM_5mm_TFQMR.mat';
file_name_hffs_e = 'Sphere_15cmrad_5mm_HFSS_Efield.fld';
file_name_hffs_h = 'Sphere_15cmrad_5mm_HFSS_Hfield.fld';

fname = fullfile(path_to_results, file_name);

%% load file

load(fname);

%% pick solution

% implicit basis-based
Sol_basis = SOL_marie2{1};


%% Post processing

axial_to_abs = @(E, slice) sqrt(E(:,:,slice,1).^2 + E(:,:,slice,2).^2 + E(:,:,slice,3).^2);
sagit_to_abs = @(E, slice) squeeze(sqrt(E(:,slice,:,1).^2 + E(:,slice,:,2).^2 + E(:,slice,:,3).^2));
coron_to_abs = @(E, slice) sqeeze(sqrt(E(slice,:,:,1).^2 + E(slice,:,:,2).^2 + E(slice,:,:,3).^2));

E_tot_basis = Sol_basis.E_tot_;

E_basis_ch1 = E_tot_basis(:,:,:,:,1);
E_basis_ch2 = E_tot_basis(:,:,:,:,2);


% compute resulting field
E_basis = E_basis_ch1 + 1i * E_basis_ch2;

E_basis = E_basis(3:63, 3:63,:,:);

idx = find(E_basis == 0);


%% read and parse hffs file 
addpath(path_to_hffs);

E_hffs = read_fld_field(file_name_hffs_e);

E_hffs = E_hffs(:,:,4:58,:);

E_hffs(idx) = 0;

%% 

slice_axial = 31;
slice_sagit = 31;
slice_coron = 31;

E_basis_abs_axial = axial_to_abs(E_basis, slice_axial);
E_basis_abs_sagit = sagit_to_abs(E_basis, slice_sagit);
E_basis_abs_coron = sagit_to_abs(E_basis, slice_coron);

figure(); imagesc(20 * log10(abs(E_basis_abs_axial)));
figure(); imagesc(20 * log10(abs(E_basis_abs_sagit)));
figure(); imagesc(20 * log10(abs(E_basis_abs_coron)));

% 
% figure(); imagesc(abs(E_basis_abs_axial));
% figure(); imagesc(abs(E_basis_abs_sagit));
% figure(); imagesc(abs(E_basis_abs_coron));
%% 

slice_axial = 30;
slice_sagit = 31;
slice_coron = 31;

E_hffs_abs_axial = axial_to_abs(E_hffs, slice_axial);
E_hffs_abs_sagit = sagit_to_abs(E_hffs, slice_sagit);
E_hffs_abs_coron = sagit_to_abs(E_hffs, slice_coron);

figure(); imagesc(20 * log10(abs(E_hffs_abs_axial)));
figure(); imagesc(20 * log10(abs(E_hffs_abs_sagit)));
figure(); imagesc(20 * log10(abs(E_hffs_abs_coron)));

% 
% figure(); imagesc(abs(E_hffs_abs_axial));
% figure(); imagesc(abs(E_hffs_abs_sagit));
% figure(); imagesc(abs(E_hffs_abs_coron));


