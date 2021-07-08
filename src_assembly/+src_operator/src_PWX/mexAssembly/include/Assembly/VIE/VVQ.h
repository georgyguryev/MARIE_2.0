#ifndef VVQ_H
#define VVQ_H

#include "EMMC.h"
#include "util.h"
#include "Gauss/gauss.h"
#include "tensors.h"
#include "VIE.h"
#include <cstring>

template<typename F>
struct GrN {
	complex<F> xx;
	complex<F> xy;
	complex<F> xz;
	complex<F> yy;
	complex<F> yz;
	complex<F> zz;
	constexpr GrN(vec3<F> const&,vec3<F> const&,reduce_type_t<F> const);
	GrN() = default;
	//#pragma omp declare simd uniform(i)
	constexpr complex<F> const& operator[](size_t i) const {
		//GrN<F> const& r = *this;
		return reinterpret_cast<complex<F> const(&)[6]>(*this)[i];
	}
	template <size_t L>
	inline void accum(F const& w,std::complex<F> (&K)[6][L],size_t const& o) const {
		K[0][o] += w*xx;
		K[1][o] += w*xy;
		K[2][o] += w*xz;
		K[3][o] += w*yy;
		K[4][o] += w*yz;
		K[5][o] += w*zz;
	}
	template <typename G>
	constexpr GrN<F>& operator+=(GrN<G> const& a) {
		xx += a.xx;
		xy += a.xy;
		xz += a.xz;
		yy += a.yy;
		yz += a.yz;
		zz += a.zz;
		return *this;
	}
	template <typename G>
	constexpr GrN<F>& operator+=(G const& a) {
		xx += a;
		xy += a;
		xz += a;
		yy += a;
		yz += a;
		zz += a;
		return *this;
	}
	template <typename G>
	constexpr GrN<F>& operator*=(GrN<G> const& a) {
		xx *= a.xx;
		xy *= a.xy;
		xz *= a.xz;
		yy *= a.yy;
		yz *= a.yz;
		zz *= a.zz;
		return *this;
	}
	template <typename G>
	constexpr GrN<F>& operator*=(G const& a) {
		xx *= a;
		xy *= a;
		xz *= a;
		yy *= a;
		yz *= a;
		zz *= a;
		return *this;
	}
};

template<typename F>
struct GrK {
	complex<F> xy;
	complex<F> xz;
	complex<F> yz;
	constexpr GrK(vec3<F> const&,vec3<F> const&,F const);
	constexpr GrK() = default;
	#pragma omp declare simd uniform(i)
	constexpr complex<F> const& operator[](size_t i) const {
		return reinterpret_cast<complex<F> const(&)[3]>(*this)[i];
	}
	template <size_t L>
	inline void accum(F const& w,std::complex<F> (&K)[3][L],size_t const& o) const {
		K[0][o] += w*xy;
		K[1][o] += w*xz;
		K[2][o] += w*yz;
	}
	template <typename G>
	constexpr GrK<F>& operator+=(GrK<G> const& a) {
		xy += a.xy;
		xz += a.xz;
		yz += a.yz;
	}
	template <typename G>
	constexpr GrK<F>& operator+=(G const& a) {
		xy += a;
		xz += a;
		yz += a;
	}
	template <typename G>
	constexpr GrK<F>& operator*=(GrK<G> const& a) {
		xy *= a.xy;
		xz *= a.xz;
		yz *= a.yz;
	}
	template <typename G>
	constexpr GrK<F>& operator*=(G const& a) {
		xy *= a;
		xz *= a;
		yz *= a;
	}
};

// e^(ix)
template<typename F>
constexpr complex<F> expmi(F x)
{
	return complex<F>(cos(x),-sin(x));
}
template <typename F,size_t N>
constexpr complexSimd<F,N> expmi(simdV<F,N> const& x) {
	complexSimd<F,N> y;
	Vc::sincos(-x,&(y.imag()),&(y.real()));
	return y;
	//return complexSimd<F,N>(cos(x),-1*sin(x));
}
//#pragma omp declare simd uniform(k0)
template<typename F>
constexpr GrN<F>::GrN(vec3<F> const& r, vec3<F> const& q,reduce_type_t<F> const k0) {
	vec3<F> const dr = r-q;
	vec3<F> const dr2 = dr*dr;
	vec3<F> const dq2(dr.x*dr.y,dr.x*dr.z,dr.y*dr.z);
	reduce_type_t<F> const k2 = k0*k0;
	F const R2 = magsq(dr);
	F const R = sqrt(R2);
	F const iR = 1/R;
	F const iR2 = iR*iR;
	F const iR3 = iR*iR2;
	F const iR4 = iR*iR3;
	F const phi = k0*R;
	//F const iR5 = iR*iR4;
	//complex<F> const Gexp = complex<F>(cos(phi),-sin(phi))*((F)M_1_4PI);
	std::complex<F> Gexp = expmi(phi);
	Gexp *= (reduce_type_t<F>)M_1_4PI;
	complex<F> const x = Gexp*complex<F>((3*dr2.x*iR2-(1+k2*dr2.x))*iR3,(3*k0*dr2.x*iR2-k0)*iR2);
	complex<F> const y = Gexp*complex<F>((3*dr2.y*iR2-(1+k2*dr2.y))*iR3,(3*k0*dr2.y*iR2-k0)*iR2);
	complex<F> const z = Gexp*complex<F>((3*dr2.z*iR2-(1+k2*dr2.z))*iR3,(3*k0*dr2.z*iR2-k0)*iR2);
	xx = -y-z;
	yy = -x-z;
	zz = -x-y;
	xy = Gexp*complex<F>((3*dq2.x*iR2-k2*dq2.x)*iR3,3*k0*dq2.x*iR4);
	xz = Gexp*complex<F>((3*dq2.y*iR2-k2*dq2.y)*iR3,3*k0*dq2.y*iR4);
	yz = Gexp*complex<F>((3*dq2.z*iR2-k2*dq2.z)*iR3,3*k0*dq2.z*iR4);
}

template<typename F>
constexpr GrK<F>::GrK(vec3<F> const& r, vec3<F> const& q, F const k0) {
	vec3<F> const dr = r-q;
	vec3<F> const dq2(dr.x*dr.y,dr.x*dr.z,dr.y*dr.z);
	F const R2 = magsq(dr);
	F const R = sqrt(R2);
	F const iR = 1/R;
	F const iR2 = iR*iR;
	F const iR3 = iR*iR2;
	F const phi = k0*R;
	//F const iR5 = iR*iR4;
	//complex<F> const G = complex<F>(cos(phi),-sin(phi))*(((F)M_1_4PI)*iR3)*complex<F>(1,phi);
	std::complex<F> G = expmi(phi);
	G *= std::complex<F>(1,phi);
	G *= ((-(reduce_type_t<F>)M_1_4PI)*iR3);
	xy = G*dr.x;
	xz = G*dr.y;
	yz = G*dr.z;
}

constexpr Size3 ind2sub3(Size3 const& D,unsigned int i) {
	unsigned int x = i % D.x;
	i /= D.x;
	unsigned int y = i % D.y;
	unsigned int z = i / D.y;
	return Size3(x,y,z);
}
constexpr unsigned int sub2ind3(Size3 const& D,Size3 const& s) {
	return s.x+D.x*(s.y+D.y*s.z);
}

#include <thread>
#include <vector>
template <typename F,bool pwl,int Nq,VIE::OpType O>
struct VVQFunctor {
	Size3 Ndom;
	Size3 Neval;
	F k0;
	vec3<F> dr;
	vec3<F> drd2;
	vec3<F> dri;
	F J;
	unsigned int N;
	unsigned int T;
	constexpr VVQFunctor(Size3 const& Ndom,Size3 const& Neval,vec3<F> const& dr,F const k0,unsigned int T) : 
	Ndom(Ndom),Neval(Neval),k0(k0),dr(dr),drd2(dr/2),dri(1/dr),J(pow(prod(drd2),2)),N(prod(Neval)),T(T) { }	
	template <int N,int s = 0,int e = ipow<6>(Nq)>
	inline std::enable_if_t<(N > 0 && s+N <= e)> eval(vec3<F> const& c,VKernel<std::complex<F>,pwl,O>& Ko) const __restrict {
		using simd_int = simdV<int,N>;
		using simd_F = simdV<F,N>;
		using GrT = std::conditional_t<VIE::isN<O>,GrN<simdV<F,N>>,GrK<simdV<F,N>>>;
		std::complex<simdV<F,N>> K[VIE::Ng<O>][VIE::Ni<pwl>] {};
		constexpr std::array<int,6> Dq{Nq,Nq,Nq,Nq,Nq,Nq};
		constexpr vec3<F> d(0,0,0);
		for(int q = s;q < e;q+=N) {
			simd_int lq = q+simd_int::IndexesFromZero();
			std::array<simd_int,6> const Q = ind2sub(lq,Dq);
		//	std::array<simd_int,6> const Q = ind2sub<simd_int,(unsigned int)Nq,6>(lq);
			simd_int const& r1 = Q[5];
			simd_int const& r2 = Q[4];
			simd_int const& r3 = Q[3];
			simd_int const& q1 = Q[2];
			simd_int const& q2 = Q[1];
			simd_int const& q3 = Q[0];
			simd_F w = Gauss<F,Nq>::w(r1,r2,r3,q1,q2,q3);
			vec3<simd_F> const r = Gauss<F,Nq>::x_v3(r1,r2,r3);
			vec3<simd_F> const p = Gauss<F,Nq>::x_v3(q1,q2,q3);
			vec3<simd_F> const rg = c+drd2*r;
			vec3<simd_F> const pg = d+drd2*p;
			GrT const g(rg,pg,k0);
			simd_F T[4] {simd_F(1),r.x/2,r.y/2,r.z/2};
			simd_F B[4] {simd_F(1),p.x/2,p.y/2,p.z/2};
			//for(auto k = 0;k < VIE::Ni<pwl>;k++) {
			//	simd_F wTBl = w*T[VIE::Il<pwl>[k]]*B[VIE::Ip<pwl>[k]];
			//	for(auto l = 0;l < VIE::Ng<O>;l++) {
			//		K[l][k] += wTBl*g[l];
			//	}
			//}
			//for(auto l = 0;l < VIE::Ng<O>;l++) {
			//for(auto k = 0;k < VIE::Ni<pwl>;k++) {
			//	simd_F wTBl = w*T[VIE::Il<pwl>[k]]*B[VIE::Ip<pwl>[k]];
			//	K[l][k] += wTBl*g[l];
			//}
			//}
			for(int k = 0;k < VIE::Ni<pwl>;k++) {
				simd_F const wTBl = w*T[VIE::Il<pwl>[k]]*B[VIE::Ip<pwl>[k]];
				g.accum(wTBl,K,k);
			}
		}
		for(int l = 0;l < VIE::Ng<O>;l++) {
		for(int k = 0;k < VIE::Ni<pwl>;k++) {
			complex<F> const lK = J*K[l][k].sum();
			if(s == 0) {
				Ko[l][k] = lK;
			}
			else {
				Ko[l][k] += lK;
			}
		}
		}
	}
	template <int N,int s = 0,int e = ipow<6>(Nq)>
	inline std::enable_if_t<(N == 0 || s+N > e)> eval(vec3<F> const&,VKernel<std::complex<F>,pwl,O>&) const { }
	inline void operator()(unsigned int const s,unsigned int const e,__restrict VKernel<F,pwl,O> Kr[],__restrict VKernel<F,pwl,O> Ki[]) const {
		constexpr vec3<F> d(0,0,0);
		constexpr size_t N = Vc::Vector<F>::size();
		constexpr size_t Tq = ipow<6>(Nq);
		constexpr size_t Nr = Tq % N;
		for(unsigned int i = s;i < e;i++) {
			Size3 const I = ind2sub3(Neval,i);
			vec3<F> const c = I*dr;
			VKernel<std::complex<F>,pwl,O> K;
			eval<N,0,Tq-Nr>(c,K);
			eval<Nr,Tq-Nr,Tq>(c,K);
			size_t const m = sub2ind3(Ndom,I);
			for(auto l = 0;l < VIE::Ng<O>;l++) {
			for(auto k = 0;k < VIE::Ni<pwl>;k++) {
				Kr[m][l][k] = std::real(K[l][k]);
				Ki[m][l][k] = std::imag(K[l][k]);
			}
			}
		}
	}
	static inline void compute(Size3 const& Ndom,Size3 const& Neval,vec3<F> const& dr,F const k0,VKernel<F,pwl,O> Kr[],VKernel<F,pwl,O> Ki[],unsigned int T=std::thread::hardware_concurrency()) {
		VVQFunctor<F,pwl,Nq,O> v(Ndom,Neval,dr,k0,T);
		std::vector<std::thread> threads;
		unsigned int B = v.N / T;
		unsigned int R = v.N % T;
		unsigned int e = 0;
		for(unsigned int t = 0;t < T;t++) {
			unsigned int s = e;
			e += B+(t<R);
			threads.push_back(std::thread(v,s,e,Kr,Ki));
		}
		for(unsigned int t = 0;t < T;t++) {
			threads[t].join();
		}
	}
};

#endif
