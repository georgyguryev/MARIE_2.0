#ifndef GAUSS_TRIG
#define GAUSS_TRIG

#include "gauss.h"
#include <array>
#include <utility>

template <typename A>
struct transform_type;
template <typename F,size_t N> struct transform_type<GaussInterval<F,N>> {
	using type = std::array<F,N>;
};
template <typename F,size_t N> struct transform_type<std::array<F,N>> {
	using type = std::array<F,N>;
};
template <typename F,size_t M,size_t N,template<typename,size_t> typename A>
struct transform_type<std::array<A<F,N>,M>> {
	using type = std::array<typename transform_type<A<F,N>>::type,M>;
};
template <typename A>
using transform_t = typename transform_type<A>::type;

/* cosine */
template <typename F,size_t ...K>
constexpr std::array<F,sizeof...(K)> cos(GaussInterval<F,sizeof...(K)> const g,std::index_sequence<K...>) {
	return std::array<F,sizeof...(K)>{(F)cos(g(K))...};
}
template <typename F,size_t Nq>
constexpr std::array<F,Nq> cos(GaussInterval<F,Nq> const g) {
	return cos<F>(g,std::make_index_sequence<Nq>{});
}
template <typename F,size_t M,size_t N,size_t ...P>
constexpr tensor<F,M,N,P...> cos(GaussTensor<F,M,N,P...> const g);
template <typename F,size_t N,size_t ...P,size_t ...K>
constexpr tensor<F,sizeof...(K),N,P...> cos(GaussTensor<F,sizeof...(K),N,P...> const g,std::index_sequence<K...>) {
	return tensor<F,sizeof...(K),N,P...>{cos<F,N,P...>(g[K])...};
}
template <typename F,size_t M,size_t N,size_t ...P>
constexpr tensor<F,M,N,P...> cos(GaussTensor<F,M,N,P...> const g) {
	return cos<F,N,P...>(g,std::make_index_sequence<M>{});
}

/* sine */
template <typename F,size_t ...K>
constexpr std::array<F,sizeof...(K)> sin(GaussInterval<F,sizeof...(K)> const g,std::index_sequence<K...>) {
	return std::array<F,sizeof...(K)>{(F)sin(g(K))...};
}
template <typename F,size_t Nq>
constexpr std::array<F,Nq> sin(GaussInterval<F,Nq> const g) {
	return sin<F>(g,std::make_index_sequence<Nq>{});
}
template <typename F,size_t M,size_t N,size_t ...P>
constexpr tensor<F,M,N,P...> sin(GaussTensor<F,M,N,P...> const g);
template <typename F,size_t N,size_t ...P,size_t ...K>
constexpr tensor<F,sizeof...(K),N,P...> sin(GaussTensor<F,sizeof...(K),N,P...> const g,std::index_sequence<K...>) {
	return tensor<F,sizeof...(K),N,P...>{sin<F,N,P...>(g[K])...};
}
template <typename F,size_t M,size_t N,size_t ...P>
constexpr tensor<F,M,N,P...> sin(GaussTensor<F,M,N,P...> const g) {
	return sin<F,N,P...>(g,std::make_index_sequence<M>{});
}

/* arctan */
template <typename F>
constexpr F atan(F x) {
	return std::atan(x);
}
template <typename F,size_t M,size_t ...N>
constexpr tensor<F,M,N...> atan(tensor<F,M,N...> const g);
template <typename F,size_t ...N,size_t ...K>
constexpr tensor<F,sizeof...(K),N...> atan(tensor<F,sizeof...(K),N...> const g,std::index_sequence<K...>) {
	return tensor<F,sizeof...(K),N...>{atan<F,N...>(g[K])...};
}
template <typename F,size_t M,size_t ...N>
constexpr tensor<F,M,N...> atan(tensor<F,M,N...> const g) {
	return atan<F,N...>(g,std::make_index_sequence<M>{});
}



#endif