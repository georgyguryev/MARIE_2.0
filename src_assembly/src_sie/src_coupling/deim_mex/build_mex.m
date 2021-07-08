% DEIM++
% Windows version takes advantage of matlab MKL library
% for Linux version it is necessary to install MKL library first  
% replace MKL_PATH with a path to your MKL directory

if ispc
    
    mex -largeArrayDims -output mexDeim  -v -O   deim_mex.cpp deim_win.cpp
      
elseif ismac
    mex -largeArrayDims -output mexDeim  -v -O   deim_mex.cpp deim_mac.cpp... 
        -llapack -lblas

elseif isunix

    
    mex -largeArrayDims -output mexDeim  -v -O   deim_mex.cpp deim_linux.cpp  ...
        COMPFLAGS='$COMPFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' CXXOPTIMFLAGS='$CXXOPTIMFLAGS -fopenmp -m64'...
        -IMKL_PATH/include ...
        -LMKL_PATH/matlab_static_libs...
        -lmkl_intel_lp64   -lmkl_gnu_thread -lmkl_core -lmkl_sequential -lmkl_core -lmkl_gnu_thread -lmkl_core -lmkl_tbb_thread -lmkl_core ... 
        -lgomp -lpthread  -ldl -lm
        
else

mex  ('-largeArrayDims','-output', 'mexDeim', '-v', '-O', '-lblas', '-llapack',...
    '-L/opt/local/lib/g95/x86_64-apple-darwin16/4.2.4/','-lf95', 'deim_mex.cpp', 'deim.cpp');

end;