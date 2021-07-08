#ifndef SSQ_LIMITS_VA_H
#define SSQ_LIMITS_VA_H

template <typename F,size_t Nq>
constexpr std::array<GaussInterval<F,Nq>,4> theta_p_va {
	GaussInterval<F,Nq>(     0, M_PI_4),
	GaussInterval<F,Nq>(     0, M_PI_4),
	GaussInterval<F,Nq>(M_PI_4, M_PI_2),
	GaussInterval<F,Nq>(M_PI_4, M_PI_2)
};
template <typename F,size_t Nq>
constexpr std::array<GaussInterval<F,Nq>,4> theta_q_va {
	GaussInterval<F,Nq>(     0, M_PI_4),
	GaussInterval<F,Nq>(M_PI_4, M_PI_2),
	GaussInterval<F,Nq>(     0, M_PI_4),
	GaussInterval<F,Nq>(M_PI_4, M_PI_2)
};

template <typename F,size_t Nq>
constexpr matrix<F,4,Nq> cos_theta_p_va = cos<F,4,Nq>(theta_p_va<F,Nq>);
template <typename F,size_t Nq>
constexpr matrix<F,4,Nq> cos_theta_q_va = cos<F,4,Nq>(theta_q_va<F,Nq>);
template <typename F,size_t Nq>
constexpr matrix<F,4,Nq> sin_theta_p_va = sin<F,4,Nq>(theta_p_va<F,Nq>);
template <typename F,size_t Nq>
constexpr matrix<F,4,Nq> sin_theta_q_va = sin<F,4,Nq>(theta_q_va<F,Nq>);


#pragma omp declare simd uniform(m)
template <typename F,bool q>
constexpr F rho_limit_va(int m,F const c,F const s) {
	F const cs[2] = {c,s};
	return 2/cs[!q ? (m/2)%2 : m%2];
}
template <typename F,size_t Nq,bool fi>
constexpr GaussInterval<F,Nq> psi_va_lim(F const psiU) {
	return fi ? GaussInterval<F,Nq>(0,psiU) : GaussInterval<F,Nq>(psiU,M_PI_2);
}

template <typename F,size_t Nq,bool fi,int m>
constexpr matrix<GaussInterval<F,Nq>,Nq,Nq> const psi_va_lim();

template <typename F,size_t Nq,bool fi>
constexpr tensor<GaussInterval<F,Nq>,4,Nq,Nq> const psi_va = {
	psi_va_lim<F,Nq,fi,0>(),
	psi_va_lim<F,Nq,fi,1>(),
	psi_va_lim<F,Nq,fi,2>(),
	psi_va_lim<F,Nq,fi,3>()
};

template <typename F,size_t Nq,bool fi>
struct trig_psi_va {
	constexpr static tensor<F,4,Nq,Nq,Nq> const c {
		cos<F,4,Nq,Nq,Nq>(psi_va<F,Nq,fi>)
	};
	constexpr static tensor<F,4,Nq,Nq,Nq> const s {
		sin<F,4,Nq,Nq,Nq>(psi_va<F,Nq,fi>)
	};
};
template <typename F,size_t Nq,bool fi>
constexpr tensor<F,4,Nq,Nq,Nq> trig_psi_va<F,Nq,fi>::c;
template <typename F,size_t Nq,bool fi>
constexpr tensor<F,4,Nq,Nq,Nq> trig_psi_va<F,Nq,fi>::s;

template <typename F,size_t Nq,bool fi,int m>
constexpr GaussInterval<F,Nq> const psi_va_lim_q1q2(int q1,int q2) {
	F cThetaP = cos_theta_p_va<F,Nq>[m][q1];
	F sThetaP = sin_theta_p_va<F,Nq>[m][q1];
	F Lp = rho_limit_va<F,0>(m,cThetaP,sThetaP);
	F cThetaQ = cos_theta_q_va<F,Nq>[m][q2];
	F sThetaQ = sin_theta_q_va<F,Nq>[m][q2];
	F Lq = rho_limit_va<F,1>(m,cThetaQ,sThetaQ);
	return psi_va_lim<F,Nq,fi>(atan(Lq/Lp));
}
template <typename F,size_t Nq,bool fi,int m,size_t ...q2>
constexpr std::array<GaussInterval<F,Nq>,Nq> const psi_va_lim_q1(int q1,std::index_sequence<q2...>) {
	static_assert(sizeof...(q2) == (size_t)Nq);
	return std::array<GaussInterval<F,Nq>,Nq>{psi_va_lim_q1q2<F,Nq,fi,m>(q1,q2)...};
}

template <typename F,size_t Nq,bool fi,int m>
constexpr std::array<GaussInterval<F,Nq>,Nq> const psi_va_lim_q1(int q1) {
	return psi_va_lim_q1<F,Nq,fi,m>(q1,std::make_index_sequence<Nq>{});
}

template <typename F,size_t Nq,bool fi,int m,size_t ...q1>
constexpr matrix<GaussInterval<F,Nq>,Nq,Nq> const psi_va_lim2(std::index_sequence<q1...>) {
	static_assert(sizeof...(q1) == Nq);
	return matrix<GaussInterval<F,Nq>,Nq,Nq>{psi_va_lim_q1<F,Nq,fi,m>(q1)...};
}
template <typename F,size_t Nq,bool fi,int m>
constexpr matrix<GaussInterval<F,Nq>,Nq,Nq> const psi_va_lim() {
	return psi_va_lim2<F,Nq,fi,m>(std::make_index_sequence<Nq>{});
}

#endif
