#ifndef DOMAIN_H
#define DOMAIN_H

#include <complex>
#include <array>
#include <initializer_list>
#include <type_traits>
#include "util.h"
using std::complex;
//using std::common_type_t;

template <typename I,size_t N,size_t ...K>
constexpr std::array<I,N-1> extract(std::array<I,N> const& i,std::index_sequence<K...>) {
	return std::array<I,N-1>{i[K+1]...};
}
template <typename I,size_t N>
constexpr std::array<I,N-1> extract(std::array<I,N> const& i) {
	return extract(i,std::make_index_sequence<N-1>{});
}

template <typename I>
constexpr I sub2ind(I const i,I const d) {
	return i;
}
template <typename I>
constexpr I sub2ind(std::array<I,1> const& i,std::array<I,1> const& d) {
	return sub2ind(i[0],d[0]);
}
template <typename I,size_t N>
constexpr I sub2ind(std::array<I,N> const& i,std::array<I,N> const& d) {
	I const& ri = i[1];
	I const& rd = d[1];
	//std::array<I,N-1> const& j = reinterpret_cast<std::array<I,N-1> const&>(ri);
	//std::array<I,N-1> const& e = reinterpret_cast<std::array<I,N-1> const&>(rd);
	return i[0]+d[0]*sub2ind(extract(i),extract(d));
}
template <typename I,size_t N,typename ...J>
constexpr std::common_type_t<I,J...> sub2ind(std::array<I,N> const& i,J const... j) {
	std::array<I,N> d = {j...};
	return sub2ind(i,d);
}

template <typename I,size_t N,size_t ...K>
constexpr std::array<I,N+1> cat_indexed(I const i,std::array<I,N> const& x,std::index_sequence<K...>) {
	return std::array<I,N+1>{i,x[K]...};
}
template <typename I,size_t N>
constexpr std::array<I,N+1> cat(I const i,std::array<I,N> const& x) {
	return cat_indexed<I,N>(i,x,std::make_index_sequence<N>{});
}
	
template <typename I,typename J>
constexpr std::array<std::common_type_t<I,J>,1> ind2sub(I const i,std::array<J,1> const& d) {
	std::array<std::common_type_t<I,J>,1> s = {i};
	return s;
}
template <typename I,typename J,size_t N>
constexpr std::array<std::common_type_t<I,J>,N> ind2sub(I const i,std::array<J,N> const& d) {
	typedef std::common_type_t<I,J> K;
	K const s0 = i % d[0];
	K const j = i / d[0];
	//K const& e0 = d[1];
	//std::array<J,N-1> const& e = reinterpret_cast<std::array<I,N-1> const&>(e0);
	std::array<K,N> s = cat<K,N-1>(s0,ind2sub(j,extract(d)));
	return s;
}
template <typename I,typename ...J,typename = std::enable_if_t<(sizeof...(J) > 1)>>
constexpr std::array<std::common_type_t<I,J...>,sizeof...(J)> ind2sub(I const i,J const... l) {
	using ctype = std::common_type_t<I,J...>;
	std::array<ctype,sizeof...(J)> const s = {ctype(l)...};
	return ind2sub(i,s);
}
template <typename I,unsigned int N,size_t D,typename = std::enable_if_t<D == 1>>
constexpr std::array<I,D> ind2sub(I const i) {
	return std::array<I,D>{i};
}
template <typename I,unsigned int N,size_t D>
constexpr std::enable_if_t<(D > 1),std::array<I,D>> ind2sub(I const i) {
	return cat<I,D-1>(i % N,ind2sub<I,N,D-1>(i / N));
}

template <unsigned int N,typename F>
constexpr F ipow(F x) {
	switch(N) {
		case 0:
			return 1;
		case 1:
			return x;
		case 2:
			return x*x;
		default: {
			if(N % 2 == 0) {
				return ipow<2,F>(ipow<N/2,F>(x));
			}
			else {
				return x*ipow<std::max(N,1u)-1,F>(x);
			}
		}
	}
}

template <typename F>
struct reduce_type {
	typedef F type;
};
template <typename F>
using reduce_type_t = typename reduce_type<F>::type;
template <typename F,size_t N>
struct reduce_type<simdV<F,N>> {
	typedef reduce_type_t<F> type;
};
template <typename F>
struct reduce_type<std::complex<F>> {
	typedef std::complex<reduce_type_t<F>> type;
};
template <typename F> struct vec3;
template <typename T,typename = std::enable_if<!is_simd<T>::value>>
constexpr vec3<T>& reduce_sum(vec3<T>& i) {
	return i;
}
template <typename T,typename = std::enable_if<is_simd<T>::value>>
constexpr vec3<reduce_type_t<T>> reduce_sum(vec3<T> const& i) {
	return vec3<reduce_type_t<T>>(i.x.sum(),i.y.sum(),i.z.sum());
}
template <typename F>
struct vec3 {
	F x;
	F y;
	F z;
	constexpr vec3(F x,F y,F z) : x(x),y(y),z(z) { }
	constexpr vec3() = default;
	constexpr vec3(F const (&r)[3]) : vec3(r[0],r[1],r[2]) { }
	inline vec3(F const r[]) : vec3(r[0],r[1],r[2]) { }
	constexpr vec3(std::array<F,3> const& r) : vec3(r[0],r[1],r[2]) { }
	F& operator[](size_t i);
	constexpr F const& operator[](size_t) const;
	constexpr std::conditional_t<is_simd<F>::value,vec3<F>,vec3<F>&> sum() {
		return reduce_sum(*this);
	}
};
//template <typename F,size_t N>
//struct vec3<simdV<F,N>> {
//	constexpr vec3<F> sum() {
//		return vec3<F>(x.sum(),y.sum(),z.sum());
//	}
//};

using Size3 = vec3<unsigned int>;
template <typename F>
struct vec2 {
	F x;
	F y;
	constexpr vec2() = default;
	constexpr vec2(F x,F y) : x(x),y(y) { }
};
template <typename F,typename G = F>
struct sym3 {
	F xx;
	G xy;
	G xz;
	F yy;
	G yz;
	F zz;
	constexpr sym3() = default;
	constexpr sym3(F xx,G xy,G xz,F yy,G yz,F zz) : xx(xx),xy(xy),xz(xz),yy(yy),yz(yz),zz(zz) { }
	constexpr std::common_type_t<F,G> const& operator[](size_t) const;
};

template <typename T>
struct is_vec3 {
	static constexpr bool value = false;
};
template <typename F>
struct is_vec3<vec3<F>> {
	static constexpr bool value = true;
};
template <typename F> inline F& vec3<F>::operator[](size_t i) {
	vec3<F>& p = *this;
	return reinterpret_cast<F(&)[3]>(p)[i];
}
template <typename F> constexpr F const& vec3<F>::operator[](size_t i) const {
	vec3<F> const& p = *this;
	return reinterpret_cast<F const(&)[3]>(p)[i];
}
template <typename F,typename G> constexpr std::common_type_t<F,G> const& sym3<F,G>::operator[](size_t i) const {
	static_assert(std::is_same<F,G>::value,"Indexing in sym3 only supported when types match");
	sym3<F,G>& p = *this;
	return reinterpret_cast<F const(&)[6]>(p)[i];
}

template <typename F> constexpr vec3<F> operator+(vec3<F> const& p) {
	return vec3<F>(+p.x,+p.y,+p.z);
}
template <typename F> constexpr vec3<F> operator-(vec3<F> const& p) {
	return vec3<F>(-p.x,-p.y,-p.z);
}
template <typename F,typename G> constexpr vec3<common_type_t<F,G>>
operator+(vec3<F> const& p,vec3<G> const& q) {
	return vec3<common_type_t<F,G>>(p.x+q.x,p.y+q.y,p.z+q.z);
}
template <typename F,typename G> constexpr vec3<common_type_t<F,G>>
operator-(vec3<F> const& p,vec3<G> const& q) {
	return vec3<common_type_t<F,G>>(p.x-q.x,p.y-q.y,p.z-q.z);
}
template <typename F,typename G> constexpr vec3<common_type_t<F,G>>
operator*(vec3<F> const& p,vec3<G> const& q) {
	return vec3<common_type_t<F,G>>(p.x*q.x,p.y*q.y,p.z*q.z);
}
// enable_if reduces ambiguity with SFINAE
template <typename F,typename G,typename = std::enable_if_t<!is_vec3<F>::value>>
constexpr vec3<common_type_t<F,G>> operator*(F const& a,vec3<G> const& q) {
	return vec3<common_type_t<F,G>>(a*q.x,a*q.y,a*q.z);
}
template <typename F,typename G,typename = std::enable_if_t<!is_vec3<G>::value>>
constexpr vec3<common_type_t<F,G>> operator*(vec3<F> const& p,G const& a) {
	return a*p;
}
template <typename F,typename G> constexpr vec3<common_type_t<F,G>>
operator/(vec3<F> const& p,G const& a) {
	return vec3<common_type_t<F,G>>(p.x/a,p.y/a,p.z/a);
}
template <typename F,typename G> constexpr vec3<common_type_t<F,G>>
operator/(F const a,vec3<G> const& p) {
	return vec3<common_type_t<F,G>>(a/p.x,a/p.y,a/p.z);
}

#include <algorithm>
template <typename F,typename G> constexpr vec3<common_type_t<F,G>> min(vec3<F> const& a,vec3<G> const& b);
template <typename ...F> constexpr common_type_t<F...> min(F const&... a) {
	using G = common_type_t<F...>;
	return std::min({((G)a)...},std::less<G>());
}
template <typename F,typename G> constexpr vec3<common_type_t<F,G>> min(vec3<F> const& a,G const& b) {
	return vec3<common_type_t<F,G>>(
		min(a.x,b),min(a.y,b),min(a.z,b)
	);
}
template <typename F,typename G> constexpr vec3<common_type_t<F,G>> min(G const& a,vec3<F> const& b) {
	return min(b,a);
}
template <typename F,typename G> constexpr vec3<common_type_t<F,G>> min(vec3<F> const& a,vec3<G> const& b) {
	return vec3<common_type_t<F,G>>(
		min(a.x,b.x),min(a.y,b.y),min(a.z,b.z)
	);
}

template <typename F> constexpr F magsq(F const& a) {
	return a*a;
}
template <typename F> constexpr F magsq(complex<F> const& a) {
	return magsq(a.real())+magsq(a.imag());
}
template <typename F> constexpr F mag(F const& a) {
	return abs(a);
}
template <typename F> constexpr F mag(complex<F> const& a) {
	return sqrt(magsq(a));
}
template <typename F> constexpr F magsq(vec3<F> const& a) {
	return magsq(a.x)+magsq(a.y)+magsq(a.z);
}
template <typename F> constexpr F mag(vec3<F> const& a) {
	return sqrt(magsq(a));
}
template <typename F> constexpr vec3<F> recip(vec3<F> const& a) {
	return 1/a;
}

template <typename F> constexpr F prod(vec3<F> const& a) {
	return a.x*a.y*a.z;
}

template <typename F,typename G> constexpr common_type_t<F,G>
dot(vec3<F> const& a,vec3<G> const& b) {
	return a.x*b.x+a.y*b.y+a.z*b.z;
}
template <typename F,typename G> constexpr common_type_t<F,G>
dot(F const (&a)[3],G const (&b)[3]) {
	return dot(reinterpret_cast<vec3<F> const&>(a),reinterpret_cast<vec3<G> const&>(b));
}
template <typename F> constexpr sym3<F>
symouter(vec3<F> const& a) {
	return sym3<F> (
		a.x*a.x,
		a.y*a.x,
		a.z*a.x,
		a.y*a.y,
		a.z*a.y,
		a.z*a.z
	);
}
template <typename F,typename G,typename H> constexpr vec3<std::common_type_t<F,G,H>>
dot(sym3<F,G> A,vec3<H> b) {
	return vec3<common_type_t<F,G,H>> (
		A.xx*b.x+A.xy*b.y+A.xz*b.z,
		A.xy*b.x+A.yy*b.y+A.yz*b.z,
		A.xz*b.x+A.yz*b.y+A.zz*b.z
	);
}
template <typename F,typename G,typename H,size_t ...K> constexpr std::array<vec3<std::common_type_t<F,G,H>>,sizeof...(K)>
dot(sym3<F,G> A,std::array<vec3<H>,sizeof...(K)> B,std::index_sequence<K...>) {
	return std::array<vec3<std::common_type_t<F,G,H>>,sizeof...(K)>{dot(A,B[K])...};
}
template <typename F,typename G,typename H,size_t N> constexpr std::array<vec3<std::common_type_t<F,G,H>>,N>
dot(sym3<F,G> A,std::array<vec3<H>,N> B) {
	return dot(A,B,std::make_index_sequence<N>{});
}
template <typename F,typename G> constexpr vec3<std::common_type_t<F,G>>
dot(std::array<vec3<F>,3> x,vec3<G> y) {
	return vec3<std::common_type_t<F,G>>(dot(x[0],y),dot(x[1],y),dot(x[2],y));
}

template <typename F,typename G>
constexpr vec3<std::common_type_t<F,G>> cross(vec3<F> const& a,vec3<G> const& b) {
	return {a.y*b.z-a.z*b.y,
			a.z*b.x-a.x*b.z,
			a.x*b.y-a.y*b.x};
}

template <typename F,typename G>
constexpr vec3<std::common_type_t<F,G>> cross(F const (&a)[3],G const (&b)[3]) {
	return cross(reinterpret_cast<vec3<F> const&>(a),reinterpret_cast<vec3<G> const&>(b));
}

template <typename F,typename U>
struct Domain {
	F* r;
	U N;
};

template <typename F,typename U>
struct Grid {
	F* x;
	F* y;
	F* z;
	U Nx;
	U Ny;
	U Nz;
};

#endif
