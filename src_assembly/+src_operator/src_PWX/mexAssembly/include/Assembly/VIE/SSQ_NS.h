#ifndef SSQ_NS_H
#define SSQ_NS_H

#include "Gauss/gauss.h"
#include "SSQ_kernel.h"

template <typename F,int Nq,bool pwl,int KType>
inline void SS_NS(complex<F> (&I)[VIE::Ni<pwl>],SS_WI<F,pwl,KType> const& W,NSKernel<F,Nq,KType> const& K) {
	F Ir[VIE::Ni<pwl>] = {0};
	F Ii[VIE::Ni<pwl>] = {0};
	#pragma omp parallel for simd collapse(4) reduction(+:Ir[:VIE::Ni<pwl>],Ii[:VIE::Ni<pwl>])
	for(int q1 = 0;q1 < Nq;q1++) {
	for(int q2 = 0;q2 < Nq;q2++) {
	for(int r1 = 0;r1 < Nq;r1++) {
	for(int r2 = 0;r2 < Nq;r2++) {
		K(Ir,Ii,bool_t<pwl>{},q1,q2,r1,r2);
	}
	}
	}
	}
	F const J[3] {K.dr.y*K.dr.z/4,K.dr.x*K.dr.z/4,K.dr.x*K.dr.y/4};
	#pragma omp simd
	for(int b = 0;b < VIE::Ni<pwl>;b++) {
		I[b] += complex<F>(Ir[b],Ii[b])*(J[K.f/2]*J[K.g/2]*W(int_t<KType>{},b));
	}
};

#endif