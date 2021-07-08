#ifndef GAUSS_H
#define GAUSS_H

#include <complex>
#include <cmath>
#include "vec3.h"
#include "gauss1D.h"
//#include "util.h"

using std::complex;
using std::real;
using std::imag;
using std::conj;

#include <cstdlib>

#include "Vc/Vc"
template <typename F,size_t Nq> struct Gauss {
	constexpr static F w(int const q) {
		return Gauss1D<F,Nq>::w[q];
	}
	constexpr static F x(int const q) {
		return Gauss1D<F,Nq>::x[q];
	}
	template <typename ...I,typename = std::enable_if_t<(sizeof...(I) > 1)>>
	constexpr static auto w(I const&... q) {
		return prod(w(q)...);
	}
	template <typename ...I>
	constexpr static std::array<F,sizeof...(I)> x(I const&... q) {
		return std::array<F,sizeof...(I)>{Gauss1D<F,Nq>::x[q]...};
	}
	template <typename I>
	constexpr static vec3<F> x_v3(I const& a,I const& b,I const& c) {
		return vec3<F>(x(a),x(b),x(c));
	}
	template <typename I,size_t N>
	inline static Vc::SimdArray<F,N> w(Vc::SimdArray<I,N> const& a) {
		//return a.apply([](I const i) { return w(i); });
		Vc::SimdArray<F,N> o;
		o.gather(Gauss1D<F,Nq>::w,a);
		return o;
	}
	template <typename I,size_t N>
	inline static Vc::SimdArray<F,N> x(Vc::SimdArray<I,N> const& a) {
		//return a.apply([](I const& i) { return x(i); });
		Vc::SimdArray<F,N> o;
		o.gather(Gauss1D<F,Nq>::x,a);
		return o;
	}
	template <typename I,size_t N>
	inline static vec3<Vc::SimdArray<F,N>> x_v3(Vc::SimdArray<I,N> const& a,Vc::SimdArray<I,N> const& b,Vc::SimdArray<I,N> const& c) {
		return vec3<Vc::SimdArray<F,N>>(x(a),x(b),x(c));
	}
};

template <typename F,size_t Nq>
struct GaussInterval {
	F a;
	F b;
	F m;
	F d;
	constexpr GaussInterval() = default;
	constexpr GaussInterval(F const& a,F const& b) : a(a),b(b),m((b+a)/2),d((b-a)/2) { }
	constexpr GaussInterval(vec2<F> const& a) : GaussInterval(a.x,a.y) { }
	constexpr GaussInterval(F const (&L)[2]) : GaussInterval(L[0],L[1]) { }
	constexpr F operator()(size_t q) const {
		return d*Gauss1D<F,Nq>::x[q]+m;
	}
};

template <typename ...F,size_t Nq>
constexpr std::common_type_t<F...> jacobian(GaussInterval<F,Nq> const&... g) {
	return prod(g.d...);
}

template <typename F,size_t Nq> struct GaussRWG {
	constexpr static size_t Nq2D = Nq*Nq;
	constexpr static vec3<F> x(int const i) {
		return x(i / Nq,i % Nq);
	}
	constexpr static vec3<F> x(int const i,int const j) {
		F const xi = Gauss<F,Nq>::x(i);
		F const xj = Gauss<F,Nq>::x(j);
		vec3<F> x { };
		x.y = (1+xj)/2;
		x.z = (1-x.y)*(1+xi)/2;
		x.x = 1-x.y-x.z;
		return x;
	}
	constexpr static F w(int const i) {
		return w(i / Nq,i % Nq);
	}
	constexpr static F w(int const i,int const j) {
		return Gauss<F,Nq>::w(i,j)*(1-Gauss<F,Nq>::x(j))/8;
	}
};

#endif
