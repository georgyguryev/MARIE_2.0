#ifndef SSQ_EA_H
#define SSQ_EA_H

#include "Gauss/gauss.h"
#include "vec3.h"
#include <array>

template <typename F>
constexpr tensor<F,2,2> simplex_ea(F const U,F const lambda,F const cTheta,F const sTheta,F const cPsi,F const sPsi) {
	return tensor<F,2,2>{
		U,lambda*sPsi-1,
		lambda*cPsi*cTheta-U,lambda*cPsi*sTheta-1
	};
}

template <typename F>
constexpr tensor<F,2,2,2> simplex_ea(F const (&U)[2],F const (&lambda)[2],F const cTheta,F const sTheta,F const cPsi,F const sPsi) {
	return tensor<F,2,2,2>{
		simplex_ea(U[0],lambda[0],cTheta,sTheta,cPsi,sPsi),
		simplex_ea(U[1],lambda[1],cTheta,sTheta,cPsi,sPsi)
	};
}

template <typename F,int Nq,bool pwl,int KType>
void SS_EA(complex<F> (&I)[VIE::Ni<pwl>],SS_WI<F,pwl,KType> const& W,VIEKernel<F,KType> const& v) {
	F Ir[VIE::Ni<pwl>] = {0};
	F Ii[VIE::Ni<pwl>] = {0};
	#pragma omp parallel for simd collapse(5) reduction(+:Ir[:VIE::Ni<pwl>],Ii[:VIE::Ni<pwl>]) default(shared)
	for(int i = 0;i < 2;i++) {
	for(int q1 = 0;q1 < Nq;q1++) {
	for(int q2 = 0;q2 < Nq;q2++) {
	for(int q3 = 0;q3 < Nq;q3++) {
	for(int q4 = 0;q4 < Nq;q4++) {
		int const m = M_ea<false>[i];
		GaussInterval<F,Nq> const thetaL = theta_ea<F,Nq>::Int[m/2];
		GaussInterval<F,Nq> const psiL = psi_ea<F,Nq>::Int[m][q1];
		GaussInterval<F,Nq> const uL = U_ea<F,Nq,0>::Int;
		GaussInterval<F,Nq> const lambdaL = lambda_ea<F,Nq,0>::Int[i][q1][q2][q3];
		F const U = uL(q3);
		F const lambda = lambdaL(q4);
		F const cPsi = psi_ea<F,Nq>::c[m][q1][q2];
		tensor<F,2,2> const Upq = simplex_ea(
			U,lambda,
			theta_ea<F,Nq>::c[m/2][q1],
			theta_ea<F,Nq>::s[m/2][q1],
			psi_ea<F,Nq>::c[m][q1][q2],
			psi_ea<F,Nq>::s[m][q1][q2]
		);
		F const w = Gauss<F,Nq>::w(q1,q2,q3,q4)*
		            thetaL.d*psiL.d*uL.d*lambdaL.d*
		            cPsi*ipow<2>(lambda);
		v.template kernel<pwl>(w,Upq,Ir,Ii);
	}
	}
	}
	}
	}
	#pragma omp parallel for simd collapse(5) reduction(+:Ir[:VIE::Ni<pwl>],Ii[:VIE::Ni<pwl>]) default(shared)
	for(int i = 0;i < 6;i++) {
	for(int q1 = 0;q1 < Nq;q1++) {
	for(int q2 = 0;q2 < Nq;q2++) {
	for(int q3 = 0;q3 < Nq;q3++) {
	for(int q4 = 0;q4 < Nq;q4++) {
		int const m = M_ea<true>[i];
		GaussInterval<F,Nq> const thetaL = theta_ea<F,Nq>::Int[m/2];
		GaussInterval<F,Nq> const psiL = psi_ea<F,Nq>::Int[m][q1];
		GaussTensor<F,2,Nq> const UL = {
			U_ea<F,Nq,1>::Int[i][q1][q2][0],
			U_ea<F,Nq,1>::Int[i][q1][q2][1]
		};
		F const U[2] = {UL[0](q3),UL[1](q3)};
		GaussTensor<F,2,Nq> const lL = {
			lambda_ea<F,Nq,1>::Int[i][q1][q2][q3][0],
			lambda_ea<F,Nq,1>::Int[i][q1][q2][q3][1]
		};
		F const L[2] = {lL[0](q4),lL[1](q4)};
		F const cPsi = psi_ea<F,Nq>::c[m][q1][q2];
		tensor<F,2,2,2> Upq = simplex_ea<F>(
			U,L,
			theta_ea<F,Nq>::c[m/2][q1],
			theta_ea<F,Nq>::s[m/2][q1],
			psi_ea<F,Nq>::c[m][q1][q2],
			psi_ea<F,Nq>::s[m][q1][q2]
		);
		F const w0 = Gauss<F,Nq>::w(q1,q2,q3,q4)*thetaL.d*psiL.d*cPsi;
		F const w[2] {
			w0*UL[0].d*lL[0].d*L[0]*L[0],
			w0*UL[1].d*lL[1].d*L[1]*L[1]
		};
		v.template kernel<pwl>(w[0],Upq[0],Ir,Ii);
		v.template kernel<pwl>(w[1],Upq[1],Ir,Ii);
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
