
% ROOT PATH == MARIE folder
ROOT_PATH = pwd;


% define paths to build mex files 
PATH_SIE_ASSEMBLY = fullfile('.','src_assembly', 'src_sie', 'src_assembly', 'src_mexdirect_ws_rwg');
PATH_C2B_ASSEMBLY = fullfile('.','src_assembly','+src_coupling', 'mex_Assembly');
PATH_PWX_ASSEMBLY = fullfile('.','src_assembly','+src_operator', 'src_PWX', 'mexAssembly', 'compile');
PATH_DIRECTFN_VIE = fullfile('.','src_assembly','+src_operator', 'src_PWX','singular');

%
BUILD_FILE_SIE_ASSEMBLY = 'mexdirect_ws_rwg_build';
BUILD_FILE_C2B_ASSEMBLY = 'build_mexAssemble_rwg_coupling';
BUILD_FILE_Linear       = 'linear_build';
BUILD_FILE_PWX_ASSEMBLY = 'compile_vie';

%%  build mex files for SIE Assembly 

cd(PATH_SIE_ASSEMBLY);
run(BUILD_FILE_SIE_ASSEMBLY);
cd(ROOT_PATH);

%% build coupling operators 

cd(PATH_C2B_ASSEMBLY);
run(BUILD_FILE_C2B_ASSEMBLY);
cd(ROOT_PATH);

%% build mex functions for PWX assembly

try
    cd(PATH_PWX_ASSEMBLY);
    run(BUILD_FILE_PWX_ASSEMBLY);
    cd(ROOT_PATH);
catch
    cd(ROOT_PATH);
    cd(PATH_DIRECTFN_VIE);
    run(BUILD_FILE_Linear)
    cd(ROOT_PATH);
end

