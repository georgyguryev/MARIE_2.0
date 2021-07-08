#ifndef SSQ_H
#define SSQ_H

#include "util.h"
#include "Gauss/gauss.h"
#include "tensors.h"
#include "points4D.h"
#include "Gauss/gaussTrig.h"
#include "SSQ_Weights.h"
#include "SSQ_kernel.h"
#include "SSQ_Limits.h"
#include "SSQ_NS.h"
#include "SSQ_VA.h"
#include "SSQ_EA.h"
#include "SSQ_ST.h"

template <typename F,bool pwl,size_t Nq,int KType>
void SSQ_Dispatch(int const m,int const f,int const g,SS_WI<F,pwl,KType> const& W,
F const k0,vec3<F> const& dr,vec3<F> const& dri,vec3<F> const& r_m,
vec3<F> const (&ord)[6][6][7],complex<F> (&I)[VIE::Ni<pwl>]) {
	if(W.Z[KType % 4])
		return;
	switch(SSAdj[m][f][g]) {
		case NSI: {
			constexpr vec3<F> const N[6] = {
				vec3<F>{-1, 0, 0},vec3<F>{+1, 0, 0},
				vec3<F>{ 0,-1, 0},vec3<F>{ 0,+1, 0},
				vec3<F>{ 0, 0,-1},vec3<F>{ 0, 0,+1}
			};
			vec3<F> const C[2] = {r_m+N[f]*dr/2,vec3<F>{0,0,0}+N[g]*dr/2};
			Gauss2DSS<F,Nq> const Gpq[2] = {
				Gauss2DSS<F,Nq>(dr,C[0]),
				Gauss2DSS<F,Nq>(dr,C[1])
			};
			NSKernel<F,Nq,KType> K(k0,dr,dri,Gpq,f,g);
			SS_NS<F,Nq,pwl,KType>(I,W,K);
		}
		break;
		case SIT: {
			VIEKernel<F,KType> const v = VIEKernelfromOrd<F,KType,SIT>(k0,dri,f,g,ord[f][g]);
			SS_ST<F,Nq,pwl,KType>(I,W,v);
		}
		break;
		case EAC:
		case EAO: {
 			VIEKernel<F,KType> const v = VIEKernelfromOrd<F,KType,EAC>(k0,dri,f,g,ord[f][g]);
			SS_EA<F,Nq,pwl,KType>(I,W,v);
		}
		break;
		case VAC:
		case VAO: {
 			VIEKernel<F,KType> const v = VIEKernelfromOrd<F,KType,VAC>(k0,dri,f,g,ord[f][g]);
			SS_VA<F,Nq,pwl,KType>(I,W,v);
		}
		break;
	}
}

template <typename F,bool pwl,size_t Nq,VIE::OpType O>
void SSQ(F const& k0,Size3 const& dims,vec3<F> const& dr,VKernel<F,pwl,O> Kr[],VKernel<F,pwl,O> Ki[]) {
	std::array<size_t,3> const D { dims.x, dims.y, dims.z };
	vec3<F>	const dri = 1/dr;
	for(size_t m = 0;m < 8;m++) {
		complex<F> I[VIE::Ng<O>][VIE::Ni<pwl>] = {0};
		std::array<size_t,3> M = ind2sub(m,2,2,2);
		if(M[0] >= dims.x || M[1] >= dims.y || M[2] >= dims.z) {
			continue;
		}
		vec3<F> const r_m = {M[0]*dr.x,M[1]*dr.y,M[2]*dr.z};
		vec3<F> ord[6][6][7];
		points_mapping<F>(m,dr,ord);
		constexpr int KOff = VIE::isN<O> ? 0 : 4;
		for(int f = 0;f < 6;f++) {
		for(int g = 0;g < 6;g++) {
			for(int i = 0;i < VIE::Ng<O>;i++) {
				SS_WO<F,pwl,O> const W(dr,k0,f,g,i);
				SSQ_Dispatch<F,pwl,Nq,0+KOff>(m,f,g,W,k0,dr,dri,r_m,ord,I[i]);
				if(pwl) {
					SSQ_Dispatch<F,pwl,Nq,1+KOff>(m,f,g,W,k0,dr,dri,r_m,ord,I[i]);
					SSQ_Dispatch<F,pwl,Nq,2+KOff>(m,f,g,W,k0,dr,dri,r_m,ord,I[i]);
					SSQ_Dispatch<F,pwl,Nq,3+KOff>(m,f,g,W,k0,dr,dri,r_m,ord,I[i]);
				}
			}
		}
		}
		int k = sub2ind(M,D);
		for(int i = 0;i < VIE::Ng<O>;i++) {
			for(int l = 0;l < VIE::Ni<pwl>;l++) {
				Kr[k][i][l] = real(I[i][l]);
				Ki[k][i][l] = imag(I[i][l]);
			}
		}
	}
}

#endif
