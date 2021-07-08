#include <iostream>
#include <cmath>
#include <complex>
#include <stdexcept>
#include "mex.h"

#ifdef _OPENMP
#include <omp.h>
#endif


using namespace std;

//------------------------------------------------------------------------ //
template <class T>
int Assemble_coupling_Lop ();
//------------------------------------------------------------------------ //

// define a Gateway function
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]) {
    
    

    cout << "Hello World! \n " << endl;
}

// EOF assembly_lop_coup.cpp
