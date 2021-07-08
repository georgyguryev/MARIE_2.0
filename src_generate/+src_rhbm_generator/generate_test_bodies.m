% generate spheres with and without implants for testing 

sigma = 1000000;

%% download a template sphere

path = './data/bodies/';

load('./data/bodies/Sphere_rad_0.05_res_0.0050.mat');

fname_sphere_1x = 'Sphere_rad_0.05_res_0.0050.mat';
fname_sphere_2x = 'Sphere_rad_0.05_res_0.0025.mat';
fname_sphere_4x = 'Sphere_rad_0.05_res_0.00125.mat';

fname_sphere_cube_1x = sprintf('Sphere_Rad_5cm_res_0.005m_with_voxelized_implant_%d.mat', sigma);
fname_sphere_cube_2x = sprintf('Sphere_Rad_5cm_res_0.0025m_with_voxelized_implant_%d.mat', sigma);
fname_sphere_cube_4x = sprintf('Sphere_Rad_5cm_res_0.00125m_with_voxelized_implant_%d.mat', sigma);


[RHBM] = src_utils.properties_to_RHBM(r, rho, epsilon_r, sigma_e);
%% refine model x2 and x4

RHBM_2x = refine_sphere(RHBM,2);
RHBM_4x = refine_sphere(RHBM,4);


%% save new RHBM models

[r, rho, epsilon_r, sigma_e] = src_utils.RHBM_to_properties(RHBM);
save(fullfile(path, fname_sphere_1x), 'r', 'rho', 'epsilon_r', 'sigma_e');
clear r rho epsilon_r sigma_e;


[r, rho, epsilon_r, sigma_e] = src_utils.RHBM_to_properties(RHBM_2x);
save(fullfile(path, fname_sphere_2x), 'r', 'rho', 'epsilon_r', 'sigma_e');
clear r rho epsilon_r sigma_e;

[r, rho, epsilon_r, sigma_e] = src_utils.RHBM_to_properties(RHBM_4x);
save(fullfile(path, fname_sphere_4x), 'r', 'rho', 'epsilon_r', 'sigma_e');
clear r rho epsilon_r sigma_e RHBM RHBM_2x RHBM_4x;

%% create a sphere with a hollow cube 1x

load(fullfile(path, fname_sphere_1x));

epsilon_center = epsilon_r(10,10,10);
sigma_center   = sigma_e(10,10,10);

% introduce conducting cube inside the sphere
epsilon_r(7:13,7:13,7:13) = 1;
sigma_e(7:13,7:13,7:13)   = sigma;

% make it hollow, fill void with material
epsilon_r(9:11,9:11,9:11) = epsilon_center;
sigma_e(9:11,9:11,9:11)   = sigma_center;

save(fullfile(path, fname_sphere_cube_1x), 'r', 'rho', 'epsilon_r', 'sigma_e');

%% create a sphere with a hollow cube 2x

load(fullfile(path, fname_sphere_2x));

epsilon_center = epsilon_r(50,50,50);
sigma_center   = sigma_e(50,50,50);

% introduce conducting cube inside the sphere
epsilon_r(48:71,48:71,48:71) = 1;
sigma_e(48:71,48:71,48:71)   = sigma;

% make it hollow, fill void with material
epsilon_r(54:65,54:65,54:65) = epsilon_center;
sigma_e(54:65,54:65,54:65)   = sigma_center;

save(fullfile(path, fname_sphere_cube_2x), 'r', 'rho', 'epsilon_r', 'sigma_e');

%% create a sphere with a hollow cube 4x

load(fullfile(path, fname_sphere_4x));

epsilon_center = epsilon_r(100,100,100);
sigma_center   = sigma_e(100,100,100);

% introduce conducting cube inside the sphere
epsilon_r(96:142,96:142,96:142) = 1;
sigma_e(96:142,96:142,96:142)   = sigma;

% make it hollow, fill void with material
epsilon_r(108:130,108:130,108:130) = epsilon_center;
sigma_e(108:130,108:130,108:130)   = sigma_center;

save(fullfile(path, fname_sphere_cube_4x), 'r', 'rho', 'epsilon_r', 'sigma_e');

clear all;

