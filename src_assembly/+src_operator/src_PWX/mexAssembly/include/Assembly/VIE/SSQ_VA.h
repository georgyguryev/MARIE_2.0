#ifndef SSQ_VA_H
#define SSQ_VA_H

#include <cmath>
#include <typeinfo>

#include "Gauss/gauss.h"
#include "vec3.h"
#include "SSQ_Limits.h"

// fi = evaluate first interval, second otherwise
template <typename F,size_t Nq,bool pwl,int KType,bool fi=true>
void inline SS_VA(complex<F> (&I)[VIE::Ni<pwl>],SS_WI<F,pwl,KType> const& W,VIEKernel<F,KType> const& v) {
	F Ir[VIE::Ni<pwl>] = {0};
	F Ii[VIE::Ni<pwl>] = {0};
	#pragma omp parallel for simd collapse(5) reduction(+:Ir[:VIE::Ni<pwl>],Ii[:VIE::Ni<pwl>])
	for(int m = 0;m < 4;m++) {
	for(size_t q1 = 0;q1 < Nq;q1++) {
	for(size_t q2 = 0;q2 < Nq;q2++) {
	for(size_t q3 = 0;q3 < Nq;q3++) {
	for(size_t q4 = 0;q4 < Nq;q4++) {
		F const cThetaP = cos_theta_p_va<F,Nq>[m][q1];
		F const sThetaP = sin_theta_p_va<F,Nq>[m][q1];
		F const cThetaQ = cos_theta_q_va<F,Nq>[m][q2];
		F const sThetaQ = sin_theta_q_va<F,Nq>[m][q2];
		F const Lp = rho_limit_va<F,0>(m,cThetaP,sThetaP);
		F const Lq = rho_limit_va<F,1>(m,cThetaQ,sThetaQ);
		F const cPsi = trig_psi_va<F,Nq,fi>::c[m][q1][q2][q3];
		F const sPsi = trig_psi_va<F,Nq,fi>::s[m][q1][q2][q3];
		GaussInterval<F,Nq> const lambdaL(0,fi?Lp/cPsi:Lq/sPsi);
		F const lambda = lambdaL(q4);
		tensor<F,2,2> const Upq = {
			lambda*cPsi*cThetaP-1,lambda*cPsi*sThetaP-1,
			lambda*sPsi*cThetaQ-1,lambda*sPsi*sThetaQ-1
		};
		
		F const w = Gauss<F,Nq>::w(q1,q2,q3,q4)*ipow<3>(lambda)*sPsi*cPsi*
		            jacobian(psi_va<F,Nq,fi>[m][q1][q2],lambdaL);
		v.template kernel<pwl>(w,Upq,Ir,Ii);
	}
	}
	}
	}
	}
	#pragma omp simd
	for(int b = 0;b < VIE::Ni<pwl>;b++) {
		I[b] += ((F)ipow<2>(M_PI_4/2))*W(int_t<KType>{},b)*complex<F>(Ir[b],Ii[b]);
	}
	if(fi) {
		SS_VA<F,Nq,pwl,KType,false>(I,W,v);
	}
}

#endif