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
    if(nrhs != 6) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
                          "Seven inputs required.");
    }
    
    // check the output parameters
    if(nlhs != 2) {
      mexErrMsgIdAndTxt( "MATLAB:convec:maxlhs",
                         "Too many output arguments.");
    }
    
    // get voxel coordinates 
    double * Scoord = mxGetDoubles(prhs[0]);
    
    // get location of sources 
    double * r_pts = mxGetDoubles(prhs[1]);
    
    // get VIE quadrature 
    double * vie_quads = mxGetDoubles(prhs[2]);
    
    // get voxel resolution
    const double res   = mxGetScalar(prhs[3]);
        
    // get voxel basis type  
    mwSize N_l  = mxGetScalar(prhs[4]);
    
    // if vie_quads has only one quadrature point - automatically switch to PWC
    if (1==vie_quads[0]) {
        N_l = 1;
    } 
    
    // get free-space wavenumber
    double ko   = mxGetScalar(prhs[5]);
    
    // get 
    mwSize N_vox  = mxGetN(prhs[0]);
    mwSize N_axes = mxGetM(prhs[0]);
    
    // get number of sources
    mwSize N_points = mxGetM(prhs[1]);
    
    // get size of 
    mwSize N_dofs = N_vox * N_axes * N_l;
    
    // allocate memory for output
    plhs[0] = mxCreateDoubleMatrix(N_dofs, N_axes * N_points, mxCOMPLEX);
    
    // allocate memory for magnetic field
    plhs[1] = mxCreateDoubleMatrix(N_dofs, N_axes * N_points, mxCOMPLEX);
    
    // get a pointer to output data
    mxComplexDouble * Z_Nop_out = mxGetComplexDoubles(plhs[0]);
    mxComplexDouble * Z_Kop_out = mxGetComplexDoubles(plhs[1]);
    
    // instantiate Factory class
    Assemble_p2v_coupling_factory factory;
    
    // call factory method that returns appropriate type for 
    Base_Assemble_p2v_coupling * coupling_assembler = factory.get_assembler(N_l);
                    
    // call appropriate assembly method
    coupling_assembler->assemble(Z_Nop_out, Z_Kop_out, Scoord, r_pts, vie_quads, res, N_l, N_vox, N_points, ko);
    
}