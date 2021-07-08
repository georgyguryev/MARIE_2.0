
%% setup compilation flags

% function compile_coupling(verbose,debug)
% 	if nargin < 1
% 		verbose = false;
% 	elseif ~islogical(verbose)
% 		verbose = false;
% 	elseif ~isscalar(verbose)
% 		verbose = verbose(1);
% 	end
% 	if nargin < 2
% 		debug = false;
% 	elseif ~islogical(debug)
% 		debug = false;
% 	elseif ~isscalar(debug)
% 		debug = debug(1);
% 	end

debug = false;
verbose = false;

FLAGS=cell(0,1);
FLAGS{end+1}='-funroll-loops -fopenmp -std=c++17 -lgomp -fPIC';
FLAGS{end+1}='-fconstexpr-depth=100000 -ftemplate-depth=100000';
FLAGS{end+1}='--param large-function-growth=2000 --param inline-unit-growth=2000 -finline-limit=3000';
FLAGS{end+1}='-fno-signed-zeros -fno-signaling-nans -fno-trapping-math -fassociative-math';
FLAGS{end+1}='-Wall -Wopenmp-simd -Wvector-operation-performance -Winline';
FLAGS{end+1}='-Wunused-const-variable';
FLAGS{end+1}='-march=native -O3 -DUSE_SIMD -DASSEMBLE_BOTH';
FLAGS=strjoin(FLAGS);
CXXFLAGS = ['CXXFLAGS="',FLAGS,'"'];
LDFLAGS = ['LDFLAGS="',FLAGS,' -t'];
LDFLAGS = [LDFLAGS,'"'];
CXXOPTIMFLAGS='CXXOPTIMFLAGS="-DNDEBUG"';
SRCDIR=mfilename('fullpath');
[SRCDIR,~,~]=fileparts(SRCDIR);
BINDIR=fullfile(SRCDIR,'bin');
LIBDIR=fullfile(SRCDIR,'lib');
INCDIR=fullfile(SRCDIR,'include');
if ~exist(BINDIR,'dir')
    mkdir(BINDIR);
end
IPATH=fullfile(SRCDIR,'include');
mexArgs = {CXXFLAGS,CXXOPTIMFLAGS,LDFLAGS,'-outdir',BINDIR,['-I',IPATH]};
if verbose
    mexArgs = [{'-v'},mexArgs];
end
if debug
    mexArgs = [{'-g'},mexArgs];
end

%% 

mex(mexArgs{:},'-R2018a', '-output', 'Assemble_rwg_coupling_matrix',...
    fullfile(SRCDIR,'src','mexAssemble_rwg_coupling_matrix.cpp'),...
    fullfile(SRCDIR,'src', 'Assemble_rwg_coupling_matrix.cpp'));

mex(mexArgs{:},'-R2018a', '-output', 'Assemble_rwg_coupling',...
    fullfile(SRCDIR,'src','mexAssemble_rwg_coupling.cpp'),...
    fullfile(SRCDIR,'src', 'Assemble_rwg_coupling.cpp'));

mex(mexArgs{:},'-R2018a', '-output', 'Assemble_p2v_coupling',...
    fullfile(SRCDIR,'src','mexAssemble_p2v_coupling.cpp'),...
    fullfile(SRCDIR,'src', 'Assemble_p2v_coupling.cpp'));
