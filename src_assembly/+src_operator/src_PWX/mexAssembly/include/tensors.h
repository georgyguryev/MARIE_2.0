#ifndef TENSORS_H
#define TENSORS_H

#include <utility>
#include <array>

template <typename F,size_t ...N> struct tensor_type;
template <typename F,size_t N> struct tensor_type<F,N> {
	using type = std::array<F,N>;
};
template <typename F,size_t M,size_t ...N> struct tensor_type<F,M,N...> {
	using type = std::array<typename tensor_type<F,N...>::type,M>;
};
template <typename F,size_t ...N>
using tensor_type_t = typename tensor_type<F,N...>::type;

template <typename F,size_t ...N> struct gtensor_type;
template <typename F,size_t N> struct gtensor_type<F,N> {
	using type = GaussInterval<F,N>;
};
template <typename F,size_t M,size_t ...N> struct gtensor_type<F,M,N...> {
	using type = std::array<typename gtensor_type<F,N...>::type,M>;
};
template <typename F,size_t ...N>
using gtensor_type_t = typename gtensor_type<F,N...>::type;

template <typename F,size_t ...N>
using tensor = tensor_type_t<F,N...>;
template <typename F,size_t ...N>
using GaussTensor = gtensor_type_t<F,N...>;
template <typename F,size_t M,size_t N>
using matrix = tensor_type_t<F,M,N>;
template <size_t N>
using index_constant = std::integral_constant<size_t,N>;

template <typename F>
constexpr F copy_tensor(F t) {
	return t;
}
template <typename F,size_t M,size_t ...N>
constexpr tensor<F,M,N...> copy_tensor(tensor<F,M,N...> const& t);
template <typename F,size_t M,size_t ...N,size_t... K>
constexpr tensor<F,M,N...> copy_tensor(tensor<F,M,N...> const& t,std::index_sequence<K...>) {
		return tensor<F,M,N...>{copy_tensor<F,N...>(t[K])...};
}
template <typename F,size_t M,size_t ...N>
constexpr tensor<F,M,N...> copy_tensor(tensor<F,M,N...> const& t) {
	return copy_tensor<F,M,N...>(t,std::make_index_sequence<M>{});
}

template <typename F,typename G>
struct TensorGen {
	template <size_t I,size_t ...J> constexpr static
	tensor<F,I,J...> gen() {
		return TensorGen<F,G>::generate(
			std::make_index_sequence<I>{},
			index_constant<J>{}...
		);
	}
	template <size_t ...I,size_t ...J,size_t K,size_t ...L> constexpr static
	tensor<F,sizeof...(J),K,L...>
	generate(std::index_sequence<J...>,index_constant<K>,index_constant<L>...) {
		return tensor<F,sizeof...(J),K,L...>{
			TensorGen<F,G>::template generate<I...,J>(
				std::make_index_sequence<K>{},index_constant<L>{}...
			)...
		};
	}
	template <size_t ...I,size_t ...J>
	constexpr static tensor<F,sizeof...(J)> generate(std::index_sequence<J...>) {
		return std::array<F,sizeof...(J)> {
			G::template gen<I...,J>()...
		};
	}
};
template <typename F,size_t Nq,typename G>
using GaussTensorGen = TensorGen<GaussInterval<F,Nq>,G>;

#endif
