#include <iostream>
#include <cmath>

#include "Assemble_coupling.h"
#include "mex.h"

// Author: Georgy Guryev, Cambridge, MA, 2019    

using std::cout;
using std::endl;

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks you out of
     the MEX-file) */
    
    // check the number of input arguments
    if(nrhs != 11) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
                          "Eleven inputs required.");
    }
    
    // check the output parameters
    if(nlhs > 2) {
      mexErrMsgIdAndTxt( "MATLAB:convec:maxlhs",
                         "Too many output arguments.");
    }
    
    // get voxel coordinates 
    double * Scoord = mxGetDoubles(prhs[0]);
    
    // get rwg vertices
    double * r_p = mxGetDoubles(prhs[1]);
    double * r_n = mxGetDoubles(prhs[2]);
    double * r_2 = mxGetDoubles(prhs[3]);
    double * r_3 = mxGetDoubles(prhs[4]);
   
    // get SIE and VIE quadratures
    double * sie_quads = mxGetDoubles(prhs[5]);
    double * vie_quads = mxGetDoubles(prhs[6]);
    
    // get voxel resolution
    const double res   = mxGetScalar(prhs[7]);
        
    // get voxel basis type  
    mwSize N_l  = mxGetScalar(prhs[8]);
    
    // if vie_quads has only one quadrature point - automatically switch to PWC
    if (1==vie_quads[0]) {
        N_l = 1;
    }
    
    // get free-space wavenumber
    double ko   = mxGetScalar(prhs[9]);
    
    int Kop_flag = mxGetScalar(prhs[10]);
    
    // get 
    mwSize N_vox  = mxGetN(prhs[0]);
    mwSize N_axes = mxGetM(prhs[0]);
    
    // get size of 
    mwSize N_dofs = N_vox * N_axes * N_l;
    mwSize N_rwg  = mxGetN(prhs[1]);
    
    // allocate memory for output
    plhs[0] = mxCreateDoubleMatrix(N_dofs, N_rwg, mxCOMPLEX);
    
    // get a pointer to output data
    mxComplexDouble * Zbc_Nop = mxGetComplexDoubles(plhs[0]);
    mxComplexDouble * Zbc_Kop = nullptr; 
        
    if (1 == Kop_flag) {
        
        plhs[1] = mxCreateDoubleMatrix(N_dofs, N_rwg, mxCOMPLEX);
        Zbc_Kop = mxGetComplexDoubles(plhs[1]);
    }
    
    // instantiate Factory class
    Assemble_rwg_coupling_matrix_factory factory;
    
    // call factory method that returns appropriate type for 
    Base_Assemble_rwg_coupling_matrix * coupling_assembler = factory.get_assembler(N_l);
            
    // call appropriate assembly method
    coupling_assembler->assemble(Zbc_Nop, Zbc_Kop, Scoord, r_p, r_n, r_2, r_3, sie_quads, vie_quads, res, N_l, N_vox, N_rwg, ko);
    
}