#ifndef VIE_H
#define VIE_H

namespace VIE {
	enum OpType { N, K };
	template <OpType O> constexpr bool isN = (O == N);
	template <OpType O> constexpr unsigned char Ng = isN<O> ?  6 : 3;
	template <bool pwl> constexpr unsigned char Nb =    pwl ?  4 : 1;
	template <bool pwl> constexpr unsigned char Ni =    pwl ? 10 : 1;
	template<bool pwl> constexpr int Il[Nb<pwl>];
	template<bool pwl> constexpr int Ip[Nb<pwl>];
	template<> constexpr int Il<true>[10] = {0,1,2,3,0,0,0,1,1,2};//{1,2,3,4,1,1,1,2,2,3};
	template<> constexpr int Ip<true>[10] = {0,1,2,3,1,2,3,2,3,3};//{1,2,3,4,2,3,4,3,4,4};
	template<> constexpr int Il<false>[1] = {0};
	template<> constexpr int Ip<false>[1] = {0};
}

// deprecated -- do not use
// namespace NOp {
// 	template<bool pwl> constexpr int Nb = pwl ? 10 : 1;
// 	template<bool pwl> constexpr int Nl = pwl ? 4 : 1;
// 	template<bool pwl> constexpr int Il[Nb<pwl>];
// 	template<bool pwl> constexpr int Ip[Nb<pwl>];
// 	template<bool pwl> constexpr int L[Nl<pwl>][Nl<pwl>];
// 	template<> constexpr int Il<true>[10] = {0,1,2,3,0,0,0,1,1,2};//{1,2,3,4,1,1,1,2,2,3};
// 	template<> constexpr int Ip<true>[10] = {0,1,2,3,1,2,3,2,3,3};//{1,2,3,4,2,3,4,3,4,4};
// 	template<> constexpr int Il<false>[1] = {0};
// 	template<> constexpr int Ip<false>[1] = {0};
// 	template<> constexpr int L<true>[4][4] {
// 		{0,4,5,6},
// 		{4,1,7,8},
// 		{5,7,2,9},
// 		{6,8,9,3}
// 	};
// 	template<> constexpr int L<false>[1][1] {
// 		{0}
// 	};
// }

template <typename F,bool pwl>
using NKernel = tensor<F,6,VIE::Ni<pwl>>;
template <typename F,bool pwl>
using KKernel = tensor<F,3,VIE::Ni<pwl>>;
template <typename F,bool pwl,VIE::OpType O>
using VKernel = std::conditional_t<VIE::isN<O>,NKernel<F,pwl>,KKernel<F,pwl>>;

#endif