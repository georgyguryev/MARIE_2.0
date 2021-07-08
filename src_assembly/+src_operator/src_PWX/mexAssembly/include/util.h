#ifndef UTIL_H
#define UTIL_H

#include <utility>
#include <complex>
#include "Vc/Vc"
#include "Vc/traits/type_traits.h"
#include "Vc/version.h"

/*********
 * types *
 *********/
template <bool b> struct bool_t { };
template <int i> struct int_t { };
template <unsigned char c> struct uchar_t { };
template <size_t s> struct sizet_t { };
template<std::size_t... I>
using uchar_sequence = std::integer_sequence<unsigned char, I...>;
template<unsigned char N>
using make_uchar_sequence = std::make_integer_sequence<unsigned char,N>;

/*******************
 * std::common_type
 *******************/
// aliases
using std::common_type_t;
template <typename F,size_t N>
using simdV = Vc::SimdArray<F,N>;
template <typename F,size_t N>
using complexSimd = std::complex<Vc::SimdArray<F,N>>;
// is_complex
template <typename T>
struct is_complex {
	static constexpr bool value = false;
};
template <typename F>
struct is_complex<std::complex<F>> {
	static constexpr bool value = true;
};
template <typename T>
struct is_complex_simd {
	static constexpr bool value = false;
};
template <typename F,size_t N>
struct is_complex_simd<complexSimd<F,N>> {
	static constexpr bool value = true;
};
template <typename T>
struct is_simd {
	static constexpr bool value = false;
};
template <typename F,size_t N>
struct is_simd<simdV<F,N>> {
	static constexpr bool value = true;
};
template <typename F,size_t N>
struct is_simd<complexSimd<F,N>> {
	static constexpr bool value = true;
};
template <typename T>
struct base_simd_type {
	typedef void type;
};
template <typename F,size_t N>
struct base_simd_type<simdV<F,N>> {
	typedef F type;
};
template <typename F,size_t N>
struct base_simd_type<complexSimd<F,N>> {
	typedef F type;
};
template <typename T>
using base_simd_type_t = typename base_simd_type<T>::type;
// common_type_if
template <typename A,typename B,typename Enable = void>
struct common_type_if { };
template <typename A,typename B,typename Enable = void>
using common_type_if_t = typename common_type_if<A,B,Enable>::type;
// enable_if
template <typename F,typename T=void>
using enable_if_arit_t = std::enable_if_t<std::is_arithmetic<F>::value,T>;
template <typename F,typename T=void>
using enable_if_complex_t = std::enable_if_t<is_complex<F>::value,T>;
template <typename F,typename T=void>
using enable_if_not_complex_t = std::enable_if_t<!is_complex<F>::value,T>;
template <typename F,typename T=void>
using enable_if_complex_simd_t = std::enable_if_t<is_complex_simd<F>::value,T>;
template <typename F,typename T=void>
using enable_if_not_complex_simd_t = std::enable_if_t<!is_complex_simd<F>::value,T>;
template <typename F,typename G,typename T=void>
using enable_if_same_t = std::enable_if_t<std::is_same<F,G>::value,T>;
template <typename F,typename G,typename T=void>
using enable_if_not_same_t = std::enable_if_t<!std::is_same<F,G>::value,T>;

/*** real-valued ***/
// common_type_if
// real SIMD, real scalar
template <typename F,typename G,size_t N>
struct common_type_if<simdV<F,N>,G,enable_if_arit_t<G>> {
	typedef simdV<std::common_type_t<F,G>,N> type;
};
// common_type
template <typename F,typename G,size_t N>
struct std::common_type<simdV<F,N>,simdV<G,N>> {
	typedef simdV<std::common_type_t<F,G>,N> type;
};
template <typename F,typename G,size_t N>
struct std::common_type<simdV<F,N>,G> {
	typedef common_type_if_t<simdV<F,N>,G,enable_if_arit_t<G>> type;
};
template <typename F,typename G,size_t N>
struct std::common_type<F,simdV<G,N>> {
	typedef std::common_type_t<simdV<G,N>,F> type;
};

/*** complex-valued ***/
// common_type_if
// real SIMD, complex scalar
template <typename F,typename G,size_t N>
struct common_type_if<simdV<F,N>,std::complex<G>,enable_if_arit_t<G>> {
	typedef std::complex<common_type_t<G,simdV<F,N>>> type;
};
// complex SIMD, complex scalar
template <typename F,typename G,size_t N>
struct common_type_if<complexSimd<F,N>,std::complex<G>,enable_if_arit_t<G>> {
	typedef std::complex<common_type_t<simdV<F,N>,G>> type;
};
// complex SIMD, real scalar
template <typename F,typename G,size_t N>
struct common_type_if<complexSimd<F,N>,G,enable_if_arit_t<G>> {
	typedef std::complex<common_type_t<simdV<F,N>,G>> type;
};
// common_type
// real SIMD, complex scalar
template <typename F,typename G,size_t N>
struct std::common_type<simdV<F,N>,std::complex<G>> {
	typedef common_type_if_t<simdV<F,N>,std::complex<G>,enable_if_arit_t<G>> type;
};
template <typename F,typename G,size_t N>
struct std::common_type<std::complex<F>,simdV<G,N>> {
	typedef common_type_t<simdV<G,N>,std::complex<F>> type;
};
// complex SIMD, real scalar
template <typename F,typename G,size_t N>
struct std::common_type<complexSimd<F,N>,G> {
	typedef common_type_if_t<complexSimd<F,N>,G,enable_if_arit_t<G>> type;
};
template <typename F,typename G,size_t N>
struct std::common_type<F,complexSimd<G,N>> {
	typedef std::common_type_t<complexSimd<G,N>,F> type;
};
// complex SIMD, complex SIMD
template <typename F,typename G,size_t N>
struct std::common_type<complexSimd<F,N>,complexSimd<G,N>> {
	typedef std::complex<common_type_t<simdV<F,N>,simdV<G,N>>> type;
};
// real SIMD, complex SIMD
template <typename F,typename G,size_t N>
struct std::common_type<simdV<F,N>,std::complex<simdV<G,N>>> {
	typedef std::complex<common_type_t<simdV<F,N>,simdV<G,N>>> type;
};
template <typename F,typename G,size_t N>
struct std::common_type<std::complex<simdV<G,N>>,simdV<F,N>> {
	typedef std::common_type_t<simdV<F,N>,std::complex<simdV<G,N>>> type;
};
// complex SIMD, complex scalar
template <typename F,typename G,size_t N>
struct std::common_type<complexSimd<F,N>,std::complex<G>> {
	typedef common_type_if_t<complexSimd<F,N>,std::complex<G>,enable_if_arit_t<G>> type;
};
template <typename F,typename G,size_t N>
struct std::common_type<std::complex<G>,complexSimd<F,N>> {
	typedef std::common_type_t<complexSimd<F,N>,std::complex<G>> type;
};

/****************
 * std::complex *
 ****************/
template <typename F,size_t N>
struct std::complex<Vc::SimdArray<F,N>> {
	using T = Vc::SimdArray<F,N>;
	using value_type = T;
	Vc::SimdArray<F,N> re;
	Vc::SimdArray<F,N> im;
	constexpr complex() = default;
	constexpr complex(T const& re,T const& im = Vc::SimdArray<F,N>()) : re(re),im(im) { }
	constexpr std::complex<F> operator[](size_t k) const {
		return std::complex<F>(re[k],im[k]);
	}
	constexpr std::complex<F> sum() const {
		return std::complex<F>(re.sum(),im.sum());
	}
	constexpr std::complex<F> prod() const {
		return std::complex<F>(re.prod(),im.prod());
	}
	inline T& real() {
		return re;
	}
	inline T& imag() {
		return im;
	}
	constexpr T const& real() const {
		return re;
	}
	constexpr T const& imag() const {
		return im;
	}
	template <typename B>
	constexpr std::complex<T>& operator+=(B const& b) {
		re += b;
		return *this;
	}
	template <typename B>
	constexpr std::complex<T>& operator+=(std::complex<B> const& b) {
		re += b.real();
		im += b.imag();
		return *this;
	}
	template <typename B>
	constexpr std::complex<T>& operator-=(B const& b) {
		re -= b;
		return *this;
	}
	template <typename B>
	constexpr std::complex<T>& operator-=(std::complex<B> const& b) {
		re -= b.real();
		im -= b.imag();
		return *this;
	}
	template <typename B>
	constexpr std::complex<T>& operator*=(B const& b) {
		re *= b;
		im *= b;
		return *this;
	}
	template <typename B>
	constexpr std::complex<T>& operator*=(std::complex<B> const& b) {
		T const pre = re;
		re = pre*b.real()-im*b.imag();
		im = pre*b.imag()+im*b.real();
		return *this;
	}
	template <typename B>
	constexpr std::complex<T>& operator/=(B const& b) {
		re /= b;
		im /= b;
		return *this;
	}
	template <typename B>
	constexpr std::complex<T>& operator/=(std::complex<B> const& b) {
		B const d = b.real()*b.real()+b.imag()*b.imag();
		T const pre = re;
		re = ( pre*b.real()+im*b.imag())/d;
		im = (-pre*b.imag()+im*b.real())/d;
		return *this;
	}
};
// operator*
template <typename F,typename G,size_t N,typename = std::enable_if_t<!std::is_same_v<complexSimd<F,N>,G>&&!std::is_same_v<simdV<F,N>,G>>>
constexpr std::common_type_t<complexSimd<F,N>,G>
operator*(complexSimd<F,N> const& a,G const& b) {
	std::common_type_t<complexSimd<F,N>,G> c = a;
	c *= b;
	return c;
}
template <typename F,typename G,size_t N,typename = std::enable_if_t<!std::is_same_v<complexSimd<F,N>,G>&&!std::is_same_v<simdV<F,N>,G>>>
//constexpr enable_if_not_complex_simd_t<G,std::common_type_t<complexSimd<F,N>,G>>
constexpr std::common_type_t<complexSimd<F,N>,G>
operator*(G const& b,std::complex<Vc::SimdArray<F,N>> const& a) {
	return a*b;
}
// operator+
template <typename F,typename G,size_t N>
constexpr enable_if_not_same_t<complexSimd<F,N>,G,std::common_type_t<complexSimd<F,N>,G>>
operator+(complexSimd<F,N> const& a,G const& b) {
	std::common_type_t<complexSimd<F,N>,G> c = a;
	c += b;
	return c;
}
template <typename F,typename G,size_t N>
constexpr enable_if_not_complex_simd_t<G,std::common_type_t<complexSimd<F,N>,G>>
operator+(G const& b,complexSimd<F,N> const& a) {
	return a+b;
}
// operator-
template <typename F,typename G,size_t N>
constexpr enable_if_not_same_t<complexSimd<F,N>,G,std::common_type_t<complexSimd<F,N>,G>>
operator-(complexSimd<F,N> const& a,G const& b) {
	std::common_type_t<complexSimd<F,N>,G> c = a;
	c -= b;
	return c;
}
template <typename F,typename G,size_t N>
constexpr enable_if_not_complex_simd_t<G,std::common_type_t<complexSimd<F,N>,G>>
operator-(G const& b,complexSimd<F,N> const& a) {
	return -(a-b);
}
/********
 * prod *
 ********/
template <typename F>
constexpr F prod(F const& x) {
	return x;
}
template <typename F,typename ...G>
constexpr std::common_type_t<F,G...> prod(F const& x,G const&... y) {
	return x*prod(y...);
}

/*******
 * sum *
 *******/
template <typename F>
constexpr F sum(F const& x) {
	return x;
}
template <typename F,typename ...G>
constexpr std::common_type_t<F,G...> sum(F const& x,G const&... y) {
	return x+sum(y...);
}
template <typename F,size_t N,size_t ...K>
constexpr F sum(std::array<F,N> const& x,std::index_sequence<K...>) {
	return sum(x[K]...);
}
template <typename F,size_t N>
constexpr F sum(std::array<F,N> const& x) {
	return sum(x,std::make_index_sequence<N>{});
}

/********
 * psum *
 ********/
template <typename F,size_t N,size_t K>
constexpr F psum(std::array<F,N> const& x,sizet_t<K>) {
	return K == 0 ? 0 : sum(x,std::make_index_sequence<K>{});
}

/*********
 *  cat  *
 *********/
template <typename ...F> constexpr std::array<std::common_type_t<F...>,sizeof...(F)>
cat(F const&... x) {
	return std::array<std::common_type_t<F...>,sizeof...(x)>{ x... };
}
template <typename F,size_t N> constexpr std::array<F,N>
cat(std::array<F,N> const& x) {
	return x;
}
template <typename F,size_t ...N,typename G> constexpr std::array<std::common_type_t<F,G>,sizeof...(N)+1>
cat(std::array<F,sizeof...(N)> const& x,std::index_sequence<N...>,G const& y) {
	return cat(x[N]...,y);
}
template <typename F,size_t N,typename G> constexpr std::array<std::common_type_t<F,G>,N>
cat(std::array<F,N> const& x, G const& y) {
	return cat(x,std::make_index_sequence<N>{},y);
}
template <typename F,typename G,size_t ...M,size_t ...N> constexpr std::array<std::common_type_t<F,G>,sizeof...(M)+sizeof...(N)>
cat(std::array<F,sizeof...(M)> const& x,std::index_sequence<M...>,std::array<G,sizeof...(N)> const& y,std::index_sequence<N...>) {
	return cat(x[M]...,y[N]...);
}
template <typename F,size_t M,typename G,size_t N> constexpr std::array<std::common_type_t<F,G>,M+N>
cat(std::array<F,M> const& x,std::array<G,N> const& y) {
	return cat(x,std::make_index_sequence<M>{},y,std::make_index_sequence<N>{});
}
//template <typename F,size_t ...N,typename ...G> constexpr std::array<std::common_type_t<F,G...>,sizeof...(G)+sizeof...(N)>
//cat(std::array<F,sizeof...(N)> x,std::index_sequence<N....>,G const&... y) {
//	return cat(x[N]...,y...);
//}
//template <typename F,size_t N,typename ...G> constexpr std::array<std::common_type_t<F,G...>,sizeof...(G)+N>
//cat(std::array<F,N> x,G const&... y) {
//	return cat(x,std::make_index_sequence<N>{},y...);
//}
//template <typename ...F,typename G,size_t N> constexpr std::array<std::common_type_t<F...,G>,sizeof...(F)+N>
//cat(F const&... x,std::array<G,N> y) {
//	return cat(y, x...);
//}
//template <typename F,size_t M,typename G,size_t ...N> constexpr std::array<std::common_type_t<F,G>,M+sizeof...(N)>
//cat(std::array<F,M> x, std::array<G,sizeof...(N)> y, std::index_sequence<N...>) {
//	return cat(x, y[N]...);
//}
//template <typename F,size_t M,typename G,size_t N> constexpr std::array<std::common_type_t<F,G>,M+N>
//cat(std::array<F,M> x,std::array<F,N> y) {
//	return cat(x,y,std::make_index_sequence<N>{});
//}
template <typename F,size_t M,typename G,size_t N,typename ...H,size_t ...P> constexpr std::array<std::common_type_t<F,G,H...>,M+N+sum(P...)>
cat(std::array<F,N> const& x, std::array<G,M> const& y, std::array<H,P> const& ...z) {
	return cat(cat(x,y),z...);
}

/*******
 * rep *
 *******/
template <typename F>
constexpr std::array<F,1> rep(F const& x,sizet_t<1>) {
	return std::array<F,1>{x};
}
template <typename F,size_t N>
constexpr std::array<F,N> rep(F const& x,sizet_t<N>) {
	return cat(x,rep(x,sizet_t<N-1>{}));
}

#endif
