#include <iostream>
#include "mex.h"
#include "matrix.h"
#include "deim.h"

using namespace std;


void mexFunction ( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
     
    if (1 != nrhs) {
         mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumInputs",
            "1 input required.");
    }
        
    // Get dimensions of matrix U
    const int m = mxGetM(prhs[0]);
    const int n = mxGetN(prhs[0]);
    
    // Get matrix U 
    const double * const U_real = mxGetPr(prhs[0]);
    const double * const U_imag = mxGetPi(prhs[0]);
    
    
    // allocate memory for output
    plhs[0] = mxCreateNumericMatrix(n,1,mxINT32_CLASS, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(m,n,mxREAL);
    
    // declare output array phi and incident matrix P
    int * phi = reinterpret_cast<int *>(mxGetPr(plhs[0]));
    double * P_in = mxGetPr(plhs[1]);
    
    // here we call the deim function
    deim(U_real, U_imag, phi, P_in, m, n);
     
//     mexPrintf("The Largest real element in first column is: %d ", phi[0] + 1);
     
    // the resulting vector PHI has to be incremented (int stores C++ indexes instead of MATLAB)
     
    // copy stuff to the left hand side pointers
    
}