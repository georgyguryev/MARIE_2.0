function compile_deim(external_mkl)

% DEIM++
% Windows version takes advantage of matlab MKL library
% for Linux version it is necessary to install MKL library first  
% replace MKL_PATH with a path to your MKL directory

if nargin < 1
    
      mex -largeArrayDims -output mexDeim  -v -O -lmwblas -lmwlapack deim_mex.cpp deim_win.cpp
      
else

    
    mex -largeArrayDims -output mexDeim  -v -O   deim_mex.cpp deim_linux.cpp  ...
        COMPFLAGS='$COMPFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' CXXOPTIMFLAGS='$CXXOPTIMFLAGS -fopenmp -m64'...
        -I MKL_PATH/include ...
        -L MKL_PATH/matlab_static_libs...
        -lmkl_intel_lp64   -lmkl_gnu_thread -lmkl_core -lmkl_sequential -lmkl_core -lmkl_gnu_thread -lmkl_core -lmkl_tbb_thread -lmkl_core ... 
        -lgomp -lpthread  -ldl -lm

end