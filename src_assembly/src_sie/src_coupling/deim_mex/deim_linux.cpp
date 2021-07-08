#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include "deim.h"
#include "mex.h"
#include <omp.h>

// define speudonims to avoid incompatibility 
#define MKL_Complex16 dcomplex

// this Linux version is developed for use with MKL library
#include "mkl.h"


using namespace std;

//-------------------------------------------------------------------------

// function compares two complex numbers (used in std::max_element())
bool abs_compare(dcomplex a, dcomplex b) {
    return (std::abs(a) < std::abs(b));
}

//------------------------------------------------------------------------------

// define DEIM++ routine
int Matvec( MKL_Complex16 *  U,  MKL_Complex16 *  elem, MKL_Complex16 *  res,
         int m,  int n) {
   
    if (nullptr == U || nullptr == elem || nullptr == res) {
         mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumInputs",
                         "Matvec(): empty input matri(-x/ces)");
        return -1;
    }
    
    // setup variable
    char trans = 'N';
    double alpha = 1.0;
    dcomplex beta  = 0.0;
    int inc = 1;
    
    // use blas routine to compute matvec
    cblas_zgemv(CblasColMajor,CblasNoTrans, m, n, &alpha, U, m, elem, inc, &beta, res, inc);

    return 0;
}

//--------------------------------------------------------------------------------

// define DEIM routine
int deim(const double * const U_real, const double * const U_imag, int * const phi,
        double * const P_in, const int m, const int n)
{
   // check if the input matrix is complex, throw error otherwise
   if (nullptr == U_real || nullptr == U_imag) {
       mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumInputs",
            "mexDEIM requires complex input matrix U!");
       return -1;
   }
   
   // check if the input matrix is complex, throw error otherwise
   if (nullptr == phi || nullptr == P_in) {
     
       mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumInputs",
            "nullptr output arguments");
       return -1;
   }
   
   // we assume that the input matrix is tall
   if(m < n) {
       mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumInputs",
            "m has to be  greater than n!");
       return -1;
   }
   
   // allocate memory 
    dcomplex * U        = new dcomplex[m * n];      // basis matrix 
    dcomplex * uL       = new dcomplex[m];          // current basis vector 
    dcomplex * mvp      = new dcomplex[m];          // result of matvec
    dcomplex * res      = new dcomplex[m];          // result from LU solver
    int * IPIV          = new int[n];               // permutation matrix
    dcomplex * elem     = new dcomplex[n];          // a vector of coefficients
    dcomplex * subU     = new dcomplex [n * n];     // submatrix with current deim points

    
    // construct complex matrix U out U_real and U_imag
    #pragma omp parallel for
    for (int ij = 0; ij < m * n; ++ij) {
        U [ij] = dcomplex(U_real[ij], U_imag[ij]);
    }
    
    // initialize first element with max row entry  
    phi[0] = std::distance(&U[0], std::max_element(&U[0], &U[0] + m, abs_compare));
  
    // variables required  for LAPACK zgesv_() and BLAS zgemv_()
    int  r = 1;
    char transp = 'N'; 
    
    for (int i = 1; i < n; ++i) {
        
        // variable for status check
        int info            = 0;

        // select next basis vector uL
        std::copy(&U[0] + i * m, &U[0] + (i + 1) * m, uL);
        
        // select row entries of current column vector
        for (int idx_iter = 0; idx_iter < i; ++idx_iter) {
            
            elem[idx_iter] = U[phi[idx_iter] + i * m];//uL[phi[idx_iter]];

            for (int k = 0; k < i; ++k) {
                subU[idx_iter + i * k] = U[phi[idx_iter] + m * k];
 
            } 
        }

        //  call matlab MKL lapack routine
        zgesv( &i, &r, subU, &i, IPIV, elem, &i, &info);
        
        // check current status of LU solver
        if (0 != info) {
            mexPrintf("U matrix is singular at step: %d !\n", info);
            return -1;
        }
              
        // Compute matrix-vector product 
        Matvec(U, elem, mvp, m, i);
                
        // Compute residual
        #pragma omp parallel for 
        for (int j = 0; j < m; ++j) {
            res[j] = uL[j] - mvp[j];
        }
        
        // find index index of max element in column vector
        phi[i] = std::distance( &res[0], std::max_element(&res[0], &res[0] + m, abs_compare));

    }
    
    // sort phi in the ascending order
    std::sort(phi, phi + n);
            
    // form incident Matrix Pin
    #pragma omp parallel for
    for (int j = 0; j < n; ++j) {
        
        // define incident matrix
        P_in[phi[j] + j * m] = 1.0;
        
        //convert phi to matlab format
        ++phi[j];
    }
    
    // deallocate memory
    delete [] U;
    delete [] mvp;
    delete [] uL;
    delete [] res;
    delete [] IPIV;         
    delete [] elem;
    delete [] subU;

    // return status
    return 0;   
}
