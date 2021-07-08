#ifndef SSQ_LIMITS_EA_H
#define SSQ_LIMITS_EA_H

/************
 ** theta ***
 ************/
template <typename F,size_t Nq>
struct theta_ea {
	constexpr static GaussTensor<F,4,Nq> Int {
		GaussInterval<F,Nq>(       0,  M_PI_4),
		GaussInterval<F,Nq>(  M_PI_4,  M_PI_2),
		GaussInterval<F,Nq>(  M_PI_2,3*M_PI_4),
		GaussInterval<F,Nq>(3*M_PI_4,  M_PI  )
	};
	constexpr static F delta = M_PI_4/2;
	constexpr static tensor<F,4,Nq> c = cos<F,4,Nq>(theta_ea<F,Nq>::Int);
	constexpr static tensor<F,4,Nq> s = sin<F,4,Nq>(theta_ea<F,Nq>::Int);
	constexpr static tensor<F,4,Nq> atcos = atan<F,4,Nq>(theta_ea<F,Nq>::c);
	constexpr static tensor<F,4,Nq> atsin = atan<F,4,Nq>(theta_ea<F,Nq>::s);
};
template <typename F,size_t Nq> constexpr GaussTensor<F,4,Nq> theta_ea<F,Nq>::Int;
template <typename F,size_t Nq> constexpr F theta_ea<F,Nq>::delta;
template <typename F,size_t Nq> constexpr tensor<F,4,Nq> theta_ea<F,Nq>::c;
template <typename F,size_t Nq> constexpr tensor<F,4,Nq> theta_ea<F,Nq>::s;
template <typename F,size_t Nq> constexpr tensor<F,4,Nq> theta_ea<F,Nq>::atcos;
template <typename F,size_t Nq> constexpr tensor<F,4,Nq> theta_ea<F,Nq>::atsin;

/************
 **** psi ***
 ************/
template <typename F,size_t Nq>
struct psi_ea {
	constexpr static GaussTensor<F,8,Nq,Nq> Int =
		GaussTensorGen<F,Nq,psi_ea<F,Nq>>::template gen<8,Nq>();
	constexpr static tensor<F,8,Nq,Nq> c = cos<F,8,Nq,Nq>(psi_ea<F,Nq>::Int);
	constexpr static tensor<F,8,Nq,Nq> s = sin<F,8,Nq,Nq>(psi_ea<F,Nq>::Int);
	template <size_t i,size_t k>
	static constexpr GaussInterval<F,Nq> gen() {
		static_assert(i < 8 && k < Nq,"Indices out of range");
		constexpr F atcos = theta_ea<F,Nq>::atcos[i/2][k];
		constexpr F atsin = theta_ea<F,Nq>::atsin[i/2][k];
		GaussInterval<F,Nq> Int[8] {
			GaussInterval<F,Nq>(     0, atcos ),
			GaussInterval<F,Nq>( atcos, M_PI_2),
			GaussInterval<F,Nq>(     0, atsin ),
			GaussInterval<F,Nq>( atsin, M_PI_2),
			GaussInterval<F,Nq>(     0, atsin ),
			GaussInterval<F,Nq>( atsin, M_PI_2),
			GaussInterval<F,Nq>(     0,-atcos ),
			GaussInterval<F,Nq>(-atcos, M_PI_2)
		};
		return Int[i];
	}
};
template <typename F,size_t Nq> constexpr GaussTensor<F,8,Nq,Nq> psi_ea<F,Nq>::Int;
template <typename F,size_t Nq> constexpr tensor<F,8,Nq,Nq> psi_ea<F,Nq>::c;
template <typename F,size_t Nq> constexpr tensor<F,8,Nq,Nq> psi_ea<F,Nq>::s;

/************
 ***** M ****
 ************/
template <bool two> constexpr tensor<int,!two ? 2 : 6> M_ea;
template <> constexpr tensor<int,2> M_ea<false> {0,6};
template <> constexpr tensor<int,6> M_ea<true> {1,2,3,4,5,7}; 

/************
 ***** U ****
 ************/
template <typename F,size_t Nq,bool two> struct U_ea;
template <typename F,size_t Nq>
struct U_ea<F,Nq,0> {
	constexpr static GaussInterval<F,Nq> Int = GaussInterval<F,Nq>(-1,1);
};
template <typename F,size_t Nq>
constexpr GaussInterval<F,Nq> U_ea<F,Nq,0>::Int;
template <typename F,size_t Nq>
struct U_ea<F,Nq,1> {
	constexpr static GaussTensor<F,6,Nq,Nq,2,Nq> Int =
		GaussTensorGen<F,Nq,U_ea<F,Nq,1>>::template gen<6,Nq,Nq,2>();
	template <size_t i,size_t k1,size_t k2,size_t k3> constexpr static GaussInterval<F,Nq>
	gen() {
		static_assert(k3 < 2,"Final param must be < 2");
		constexpr int m = M_ea<true>[i];
		F U = 0;
		F Vtp = theta_ea<F,Nq>::c[m/2][k1]*psi_ea<F,Nq>::c[m][k1][k2]/psi_ea<F,Nq>::s[m][k1][k2];
		F Vth = theta_ea<F,Nq>::c[m/2][k1]/theta_ea<F,Nq>::s[m/2][k1];
		switch(m) {
			case 1:
			case 3:
				U = 2*Vtp-1;
				break;
			case 2:
				U = 2*Vth-1;
				break;
			case 4:
				U = 2*Vth+1;
				break;
			case 5:
			case 7:
				U = 2*Vtp+1;
				break;
		}
		// note: flipped case 7
		return k3 == 0 ? GaussInterval<F,Nq>(-1,U) : GaussInterval<F,Nq>(U,1);
	}
};
template <typename F,size_t Nq> constexpr GaussTensor<F,6,Nq,Nq,2,Nq> U_ea<F,Nq,1>::Int;

/************
 ** lambda **
 ************/
template <typename F,size_t Nq,bool two>
struct lambda_ea;
template <typename F,size_t Nq>
struct lambda_ea<F,Nq,0> {
	constexpr static GaussTensor<F,2,Nq,Nq,Nq,Nq> Int =
		TensorGen<GaussInterval<F,Nq>,lambda_ea<F,Nq,0>>::template gen<2,Nq,Nq,Nq>();
	template <size_t i,size_t q1,size_t q2,size_t q3>
	static constexpr GaussInterval<F,Nq> gen() {
		constexpr size_t m = M_ea<0>[i];
		constexpr F cTheta = theta_ea<F,Nq>::c[m/2][q1];
		constexpr F cPsi = psi_ea<F,Nq>::c[m][q1][q2];
		constexpr F U = U_ea<F,Nq,0>::Int(q3);
		switch(m) {
			case 0:
				return GaussInterval<F,Nq>(0,(U+1)/(cTheta*cPsi));
			case 6:
				return GaussInterval<F,Nq>(0,(U-1)/(cTheta*cPsi));
		}
	}
};
template <typename F,size_t Nq> constexpr GaussTensor<F,2,Nq,Nq,Nq,Nq> lambda_ea<F,Nq,0>::Int;
template <typename F,size_t Nq>
struct lambda_ea<F,Nq,1> {
	constexpr static GaussTensor<F,6,Nq,Nq,Nq,2,Nq> Int =
		TensorGen<GaussInterval<F,Nq>,lambda_ea<F,Nq,1>>::template gen<6,Nq,Nq,Nq,2>();
	template <size_t i,size_t q1,size_t q2,size_t q3,size_t k>
	static constexpr GaussInterval<F,Nq> gen() {
		constexpr size_t m = M_ea<true>[i];
		constexpr bool f = k == 0;
		constexpr F cTheta = theta_ea<F,Nq>::c[m/2][q1];
		constexpr F sTheta = theta_ea<F,Nq>::s[m/2][q1];
		constexpr F cPsi = psi_ea<F,Nq>::c[m][q1][q2];
		constexpr F sPsi = psi_ea<F,Nq>::s[m][q1][q2];
		constexpr F V[2] {
			(U_ea<F,Nq,1>::Int[i][q1][q2][0](q3)+1)/(cTheta*cPsi),
			(U_ea<F,Nq,1>::Int[i][q1][q2][1](q3)-1)/(cTheta*cPsi)
		};
		constexpr F W[2] {
			2/sPsi,
			2/(sTheta*cPsi)
		};
		switch(m) {
			case 1:
			case 3:
				return GaussInterval<F,Nq>(0,f ? V[0] : W[0]);
			case 2:
				return GaussInterval<F,Nq>(0,f ? V[0] : W[1]);
			case 4:
				return GaussInterval<F,Nq>(0,f ? W[1] : V[1]);
			case 5:
			case 7:
				return GaussInterval<F,Nq>(0,f ? W[0] : V[1]);
		}
	}
};
template <typename F,size_t Nq> constexpr GaussTensor<F,6,Nq,Nq,Nq,2,Nq> lambda_ea<F,Nq,1>::Int;

#endif
