function compile_vie(verbose,debug)
	if nargin < 1
		verbose = false;
	elseif ~islogical(verbose)
		verbose = false;
	elseif ~isscalar(verbose)
		verbose = verbose(1);
	end
	if nargin < 2
		debug = false;
	elseif ~islogical(debug)
		debug = false;
	elseif ~isscalar(debug)
		debug = debug(1);
	end
	FLAGS=cell(0,1);
	FLAGS{end+1}='-funroll-loops -fopenmp -std=c++17 -lgomp';
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
	SRCDIR=strsplit(SRCDIR,filesep);
	SRCDIR=strjoin(SRCDIR(1:end-1),filesep);
	BINDIR=fullfile(SRCDIR,'bin');
	LIBDIR=fullfile(SRCDIR,'lib');
	INCDIR=fullfile(SRCDIR,'include');
	if ~exist(fullfile(LIBDIR,'libVc.a'),'file') || isempty(dir(fullfile(INCDIR,'Vc')))
		system(fullfile(INCDIR,'external','build.sh'));
	end
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
 	clear vie_assembly
 	mex(mexArgs{:},fullfile(SRCDIR,'src','vie_assembly.cpp'),['-L',LIBDIR],'-lVc')
end
