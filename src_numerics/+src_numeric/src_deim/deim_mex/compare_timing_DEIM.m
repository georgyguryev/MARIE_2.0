% script for timing of matlab and cpp DEIM implementations


%% define basis filename 

basis_fname = 'Basis_298MHz_Triangular_medium_sphere_1e-8.mat';
path        = fullfile('.','data', 'basis');

basis_fname = fullfile(path,basis_fname);

%% Load basis 

load (basis_fname);

%% run  DEIM

tic;

[phi_1, ~] = src_numeric.deim(U);

t_DEIM = toc;


%% run DEIM++

tic;

[Pin] = mexDeim(U);

t_DEIM_cpp = toc;

%%

fid = 1;

fprintf(fid, '\n Time spent to compute DEIM is: %f s \n',t_DEIM);
fprintf(fid, '\n Time spent to compute DEIM++ is: %f s \n',t_DEIM_cpp);

fprintf(fid, '\n The time ratio is: %f\n', t_DEIM / t_DEIM_cpp);



