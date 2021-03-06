
#include <iostream>
//#include "directfn_inline.h"
#include "classes.h"
#include <complex>
#include <math.h>
#include "directfn.h"
using namespace std;

void create_VA(const double r1[],const double r2[],const double r3[],const double r4[],const double r5[],const double r6[], const double r7[],const int N1,const int N2, const int N3, const int N4, double k0, double dx, double rq_c[],double rp_c[],double nq[], double np[], int ker_type, int l, int lp, complex<double> I[])
{
	Geometry geom;
    geom.VA(r1,r2,r3,r4,r5,r6,r7);
    geom.set_wavenumber(k0);
    geom.set_delta(dx);
	geom.set_centers(rq_c,rp_c);
	geom.set_normales(nq,np);
	geom.set_kerneltype(ker_type);
	geom.set_lp(l);
	geom.set_lq(lp);

	I[0] = quadric_ws_va(N1, N2, N3, N4, geom);
	
}


