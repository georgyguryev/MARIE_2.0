#ifndef SSQ_LIMITS_ST_H
#define SSQ_LIMITS_ST_H

template <typename F,size_t Nq>
struct psi_st {
	static constexpr GaussTensor<F,6,Nq> Int {
		GaussInterval<F,Nq>(       0,  M_PI_4),
		GaussInterval<F,Nq>(  M_PI_4,  M_PI_2),
		GaussInterval<F,Nq>(  M_PI_4,  M_PI_2),
		GaussInterval<F,Nq>(  M_PI_2,3*M_PI_4),
		GaussInterval<F,Nq>(  M_PI_2,3*M_PI_4),
		GaussInterval<F,Nq>(3*M_PI_4,  M_PI  )
	};
	static constexpr tensor<F,6,Nq> c = cos<F,6,Nq>(psi_st<F,Nq>::Int);
	static constexpr tensor<F,6,Nq> s = sin<F,6,Nq>(psi_st<F,Nq>::Int);
};
template <typename F,size_t Nq>
constexpr GaussTensor<F,6,Nq> psi_st<F,Nq>::Int;
template <typename F,size_t Nq>
constexpr tensor<F,6,Nq> psi_st<F,Nq>::c;
template <typename F,size_t Nq>
constexpr tensor<F,6,Nq> psi_st<F,Nq>::s;

template <typename F,size_t Nq>
struct U_st {
	static constexpr GaussTensor<F,6,Nq,Nq> Int = 
		GaussTensorGen<F,Nq,U_st<F,Nq>>::template gen<6,Nq>();
	template <size_t m,size_t q1>
	static constexpr GaussInterval<F,Nq> gen() {
		F cot = psi_st<F,Nq>::c[m][q1]/psi_st<F,Nq>::s[m][q1];
		switch(m) {
			case 1:
				return GaussInterval<F,Nq>(-2*cot+1,+1);
			case 2:
				return GaussInterval<F,Nq>(-1,-2*cot+1);
			case 3:
				return GaussInterval<F,Nq>(-2*cot-1,+1);
			case 4:
				return GaussInterval<F,Nq>(-1,-2*cot-1);
			default:
				return GaussInterval<F,Nq>(-1,1);
		}
	}
};
template <typename F,size_t Nq>
constexpr GaussTensor<F,6,Nq,Nq> U_st<F,Nq>::Int;

template <typename F,size_t Nq>
struct lambda_st {
	static constexpr GaussTensor<F,6,Nq,Nq,Nq> Int = 
		GaussTensorGen<F,Nq,lambda_st<F,Nq>>::template gen<6,Nq,Nq>();
	template <size_t m,size_t q1,size_t q2>
	static constexpr GaussInterval<F,Nq> gen() {
		switch(m) {
			case 0:
			case 1:
				return GaussInterval<F,Nq>(0,(1-U_st<F,Nq>::Int[m][q1](q2))/psi_st<F,Nq>::c[m][q1]);
			case 2:
			case 3:
				return GaussInterval<F,Nq>(0,2/psi_st<F,Nq>::s[m][q1]);
			case 4:
			case 5:
				return GaussInterval<F,Nq>(0,(1+U_st<F,Nq>::Int[m][q1](q2))/-psi_st<F,Nq>::c[m][q1]);
		}
	}
};
template <typename F,size_t Nq>
constexpr GaussTensor<F,6,Nq,Nq,Nq> lambda_st<F,Nq>::Int;

template <typename F,size_t Nq>
struct rho_st {
	static constexpr GaussInterval<F,Nq> Int(F const lambda) {
		return GaussInterval<F,Nq>(0,lambda);
	}
};

template <typename F>
constexpr tensor<F,2> subtriangles_st(int const i,F const U,F const V) {
	tensor<F,4,2> const t{
		tensor<F,2>{+U,+V},
		tensor<F,2>{-V,+U},
		tensor<F,2>{-U,-V},
		tensor<F,2>{+V,-U}
	};
	return t[i];
}

template <typename F>
constexpr tensor<F,2,2> subtriangles_st(int const i,F const Up,F const Vp,F const Uq,F const Vq) {
	return tensor<F,2,2>{
		subtriangles_st<F>(i,Up,Vp),
		subtriangles_st<F>(i,Uq,Vq)
	};
}

#endif
