#include <complex>

using namespace std;


// declaire aliases
using dcomplex = std::complex <double>;
using fcomplex = std::complex <float>;

// DEIM++ declaration
int deim(const double * const U_real, const double * const U_imag, int * const phi,
        double * const P_in, const int m, const int n);

// function to compare two complex numbers by absolute value
bool abs_compare(dcomplex a, dcomplex b);

// Wrapper function around MKL zgemv() routine to compute matrix-vactor product
int Matvec( dcomplex * const U,  dcomplex * const elem, dcomplex * const res,
           int m,  int n);
