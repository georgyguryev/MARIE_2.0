#ifndef SSQ_ST_H
#define SSQ_ST_H

#include "Gauss/gauss.h"
#include "vec3.h"
#include "SSQ_kernel.h"

template <typename F,int Nq,bool pwl,int KType>
inline void SS_ST(complex<F> (&I)[VIE::Ni<pwl>],SS_WI<F,pwl,KType> const& W,VIEKernel<F,KType> const& v) {
 	F Ir[VIE::Ni<pwl>] = {0};
	F Ii[VIE::Ni<pwl>] = {0};
	#pragma omp parallel for simd collapse(6) reduction(+:Ir[:VIE::Ni<pwl>],Ii[:VIE::Ni<pwl>])
	for(int k = 0;k < 4;k++) {
	for(int m = 0;m < 6;m++) {
	for(int q1 = 0;q1 < Nq;q1++) {
	for(int q2 = 0;q2 < Nq;q2++) {
	for(int q3 = 0;q3 < Nq;q3++) {
	for(int q4 = 0;q4 < Nq;q4++) {
		GaussInterval<F,Nq> const psiL = psi_st<F,Nq>::Int[m];
		GaussInterval<F,Nq> const uL = U_st<F,Nq>::Int[m][q1];
		GaussInterval<F,Nq> const lL = lambda_st<F,Nq>::Int[m][q1][q2];
		F const lambda = lL(q3);
		GaussInterval<F,Nq> const rhoL = rho_st<F,Nq>::Int(lambda);
		F const rho = rhoL(q4);

		F const sPsi = psi_st<F,Nq>::s[m][q1];
		F const cPsi = psi_st<F,Nq>::c[m][q1];
		F const U = uL(q2);

		F const Up0 = U;
		F const Vp0 = lambda*sPsi-1;
		F const Uq0 =  rho*cPsi+Up0;
		F const Vq0 = -rho*sPsi+Vp0;
		tensor<F,2,2> const Upq = subtriangles_st(k,Up0,Vp0,Uq0,Vq0);
		
		F const w = Gauss<F,Nq>::w(q1,q2,q3,q4)*uL.d*lL.d*psiL.d*rhoL.d*sPsi*rho;
		v.template kernel<pwl>(w,Upq,Ir,Ii);
	}
	}
	}
	}
	}
	}
	#pragma omp simd
	for(int b = 0;b < VIE::Ni<pwl>;b++) {
		I[b] += complex<F>(Ir[b],Ii[b])*W(int_t<KType>{},b);
	}
}

#endif
