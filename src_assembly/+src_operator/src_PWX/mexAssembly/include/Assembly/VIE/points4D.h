#ifndef POINTS_4D_H
#define POINTS_4D_H

#include <algorithm>
#include <typeinfo>

enum SType {
	NSI,VAO,VAC,EAO,EAC,SIT
};

constexpr const char* cchar(SType s) {
	switch(s) {
		case NSI:
			return "NSI";
		case VAO:
			return "VAO";
		case VAC:
			return "VAC";
		case EAO:
			return "EAO";
		case EAC:
			return "EAC";
		case SIT:
			return "SIT";
		default:
			return "";
	}
}

constexpr SType SSAdj[8][6][6] {
	{ // 0 [1,1,1]
		{SIT,NSI,EAO,EAO,EAO,EAO},
		{NSI,SIT,EAO,EAO,EAO,EAO},
		{EAO,EAO,SIT,NSI,EAO,EAO},
		{EAO,EAO,NSI,SIT,EAO,EAO},
		{EAO,EAO,EAO,EAO,SIT,NSI},
		{EAO,EAO,EAO,EAO,NSI,SIT}
	},
	{ // 1 [2,1,1]
		{NSI,SIT,EAO,EAO,EAO,EAO},
		{NSI,NSI,NSI,NSI,NSI,NSI},
		{NSI,EAO,EAC,NSI,VAO,VAO},
		{NSI,EAO,NSI,EAC,VAO,VAO},
		{NSI,EAO,VAO,VAO,EAC,NSI},
		{NSI,EAO,VAO,VAO,NSI,EAC}
	},
	{ // 2 [1,2,1]
		{EAC,NSI,NSI,EAO,VAO,VAO},
		{NSI,EAC,NSI,EAO,VAO,VAO},
		{EAO,EAO,NSI,SIT,EAO,EAO},
		{NSI,NSI,NSI,NSI,NSI,NSI},
		{VAO,VAO,NSI,EAO,EAC,NSI},
		{VAO,VAO,NSI,EAO,NSI,EAC}
	},
	{ // 3 [2,2,1]
		{NSI,EAC,NSI,EAO,VAO,VAO},
		{NSI,NSI,NSI,NSI,NSI,NSI},
		{NSI,EAO,NSI,EAC,VAO,VAO},
		{NSI,NSI,NSI,NSI,NSI,NSI},
		{NSI,VAO,NSI,VAO,VAC,NSI},
		{NSI,VAO,NSI,VAO,NSI,VAC}
	},
	{ // 4 [1,1,2]
		{EAC,NSI,VAO,VAO,NSI,EAO},
		{NSI,EAC,VAO,VAO,NSI,EAO},
		{VAO,VAO,EAC,NSI,NSI,EAO},
		{VAO,VAO,NSI,EAC,NSI,EAO},
		{EAO,EAO,EAO,EAO,NSI,SIT},
		{NSI,NSI,NSI,NSI,NSI,NSI}
	},
	{ // 5 [2,1,2]
		{NSI,EAC,VAO,VAO,NSI,EAO},
		{NSI,NSI,NSI,NSI,NSI,NSI},
		{NSI,VAO,VAC,NSI,NSI,VAO},
		{NSI,VAO,NSI,VAC,NSI,VAO},
		{NSI,EAO,VAO,VAO,NSI,EAC},
		{NSI,NSI,NSI,NSI,NSI,NSI}
	},
	{ // 6 [1,2,2]
		{VAC,NSI,NSI,VAO,NSI,VAO},
		{NSI,VAC,NSI,VAO,NSI,VAO},
		{VAO,VAO,NSI,EAC,NSI,EAO},
		{NSI,NSI,NSI,NSI,NSI,NSI},
		{VAO,VAO,NSI,EAO,NSI,EAC},
		{NSI,NSI,NSI,NSI,NSI,NSI}
	},
	{ // 7 [2,2,2]
		{NSI,VAC,NSI,VAO,NSI,VAO},
		{NSI,NSI,NSI,NSI,NSI,NSI},
		{NSI,VAO,NSI,VAC,NSI,VAO},
		{NSI,NSI,NSI,NSI,NSI,NSI},
		{NSI,VAO,NSI,VAO,NSI,VAC},
		{NSI,NSI,NSI,NSI,NSI,NSI}
	}
};

template <typename F,int Nq>
inline vec3<F> const& points4D_const(int const& M,vec3<F> const& c,vec3<F> const& dr,int(&Iq)[2]) {
	return
	{c.x + (1-(M == 0)) * dr.x/2*Gauss1D<F,Nq>::x[Iq[0]],
	 c.y + (1-(M == 1)) * dr.y/2*Gauss1D<F,Nq>::x[Iq[M != 0]],
	 c.z + (1-(M == 2)) * dr.z/2*Gauss1D<F,Nq>::x[Iq[1]]};
}

template <typename F>
constexpr void  assign_ord(vec3<F> (&o)[7],std::initializer_list<vec3<F>>& i) {
	std::copy(std::begin(i),std::end(i),std::begin(o));
}
template <typename F,typename ...G>
constexpr void assign_ord(vec3<F> (&o)[7],vec3<G> const... I) {
	static_assert(std::is_same<std::common_type_t<G...>,F>::value);
	std::initializer_list<vec3<F>> Il = {I...};
	assign_ord(o,Il);
}
template <typename F>
using OrdV3 = std::array<vec3<F>,7>;
template <typename F>
constexpr void assign_ord(vec3<F> (&o)[7],OrdV3<F> const& i) {
	for(int k = 0;k < 7;k++) {
		o[k] = i[k];
	}
}


template <typename F,int m>
void points_mapping(vec3<F> const& dr,OrdV3<F>(&ord)[6][6]);

template <typename F>
void points_mapping(int const& m,vec3<F> const& dr,vec3<F>(&ord)[6][6][7]) {
	OrdV3<F> (&orda)[6][6] = reinterpret_cast<OrdV3<F>(&)[6][6]>(ord);
	switch(m) {
		case 0:
			points_mapping<F,0>(dr,orda);
			break;
		case 1:
			points_mapping<F,1>(dr,orda);
			break;
		case 2:
			points_mapping<F,2>(dr,orda);
			break;
		case 3:
			points_mapping<F,3>(dr,orda);
			break;
		case 4:
			points_mapping<F,4>(dr,orda);
			break;
		case 5:
			points_mapping<F,5>(dr,orda);
			break;
		case 6:
			points_mapping<F,6>(dr,orda);
			break;
		case 7:
		default:
			points_mapping<F,7>(dr,orda);
			break;
	}
}

template <typename F,int k>
constexpr vec3<F> Pord(vec3<F> const& dp) {
	switch(k) {
		case 0:
			return vec3<F>({-dp.x,-dp.y,-dp.z});
		case 1:
			return vec3<F>({+dp.x,-dp.y,-dp.z});
		case 2:
			return vec3<F>({+dp.x,+dp.y,-dp.z});
		case 3:
			return vec3<F>({-dp.x,+dp.y,-dp.z});
		case 4:
			return vec3<F>({-dp.x,-dp.y,+dp.z});
		case 5:
			return vec3<F>({+dp.x,-dp.y,+dp.z});
		case 6:
			return vec3<F>({+dp.x,+dp.y,+dp.z});
		case 7:
			return vec3<F>({-dp.x,+dp.y,+dp.z});
	}
}

template <typename F,int m>
void points_mapping(vec3<F> const& dr,OrdV3<F> (&ord)[6][6]) {
	constexpr vec3<int> M = {m % 2,(m % 4) / 2,(m / 4) % 2};
	vec3<F> const dp = dr/2;
	vec3<F> const c = M*dr;
	vec3<F> const p[8] {
		{-dp.x,-dp.y,-dp.z},
		{+dp.x,-dp.y,-dp.z},
		{+dp.x,+dp.y,-dp.z},
		{-dp.x,+dp.y,-dp.z},
		{-dp.x,-dp.y,+dp.z},
		{+dp.x,-dp.y,+dp.z},
		{+dp.x,+dp.y,+dp.z},
		{-dp.x,+dp.y,+dp.z}
	};
// 	vec3<F> const r[8] {
// 		{c+p[0]},
// 		{c+p[1]},
// 		{c+p[2]},
// 		{c+p[3]},
// 		{c+p[4]},
// 		{c+p[5]},
// 		{c+p[6]},
// 		{c+p[7]}
// 	};
	constexpr vec3<F> zvec = {0,0,0};
	switch(m) {
		case 0: // [1,1,1]
		ord[0][0] = OrdV3<F>{p[0],p[3],p[7],p[4],zvec,zvec,zvec};
		ord[0][2] = OrdV3<F>{p[1],p[5],p[4],p[0],p[3],p[7],zvec};
		ord[0][3] = OrdV3<F>{p[2],p[6],p[7],p[3],p[0],p[4],zvec};
		ord[0][4] = OrdV3<F>{p[1],p[2],p[3],p[0],p[4],p[7],zvec};
		ord[0][5] = OrdV3<F>{p[5],p[6],p[7],p[4],p[0],p[3],zvec};

		ord[1][1] = OrdV3<F>{p[1],p[2],p[6],p[5],zvec,zvec,zvec};
		ord[1][2] = OrdV3<F>{p[0],p[4],p[5],p[1],p[2],p[6],zvec};
		ord[1][3] = OrdV3<F>{p[3],p[7],p[6],p[2],p[1],p[5],zvec};
		ord[1][4] = OrdV3<F>{p[0],p[3],p[2],p[1],p[5],p[6],zvec};
		ord[1][5] = OrdV3<F>{p[7],p[4],p[5],p[6],p[2],p[1],zvec};

		ord[2][0] = OrdV3<F>{p[3],p[7],p[4],p[0],p[1],p[5],zvec};
		ord[2][1] = OrdV3<F>{p[2],p[6],p[5],p[1],p[0],p[4],zvec};
		ord[2][2] = OrdV3<F>{p[0],p[1],p[5],p[4],zvec,zvec,zvec};
		ord[2][4] = OrdV3<F>{p[2],p[3],p[0],p[1],p[5],p[4],zvec};
		ord[2][5] = OrdV3<F>{p[6],p[7],p[4],p[5],p[1],p[0],zvec};

		ord[3][0] = OrdV3<F>{p[0],p[4],p[7],p[3],p[2],p[6],zvec};
		ord[3][1] = OrdV3<F>{p[1],p[5],p[6],p[2],p[3],p[7],zvec};
		ord[3][3] = OrdV3<F>{p[3],p[2],p[6],p[7],zvec,zvec,zvec};
		ord[3][4] = OrdV3<F>{p[0],p[1],p[2],p[3],p[7],p[6],zvec};
		ord[3][5] = OrdV3<F>{p[4],p[5],p[6],p[7],p[3],p[2],zvec};

		ord[4][0] = OrdV3<F>{p[7],p[4],p[0],p[3],p[2],p[1],zvec};
		ord[4][1] = OrdV3<F>{p[5],p[6],p[2],p[1],p[0],p[3],zvec};
		ord[4][2] = OrdV3<F>{p[4],p[5],p[1],p[0],p[3],p[2],zvec};
		ord[4][3] = OrdV3<F>{p[6],p[7],p[3],p[2],p[1],p[0],zvec};
		ord[4][4] = OrdV3<F>{p[0],p[1],p[2],p[3],zvec,zvec,zvec};

		ord[5][0] = OrdV3<F>{p[0],p[3],p[7],p[4],p[5],p[6],zvec};
		ord[5][1] = OrdV3<F>{p[1],p[2],p[6],p[5],p[4],p[7],zvec};
		ord[5][2] = OrdV3<F>{p[0],p[1],p[5],p[4],p[7],p[6],zvec};
		ord[5][3] = OrdV3<F>{p[3],p[2],p[6],p[7],p[4],p[5],zvec};
		ord[5][5] = OrdV3<F>{p[4],p[5],p[6],p[7],zvec,zvec,zvec};
		break;
		case 1: // 2,1,1
		ord[0][1] = OrdV3<F>{p[1],p[2],p[6],p[5],zvec,zvec,zvec};
		ord[0][2] = OrdV3<F>{p[0],p[4],p[5],p[1],p[2],p[6],zvec};
		ord[0][3] = OrdV3<F>{p[3],p[7],p[6],p[2],p[1],p[5],zvec};
		ord[0][4] = OrdV3<F>{p[0],p[3],p[2],p[1],p[5],p[6],zvec};
		ord[0][5] = OrdV3<F>{p[4],p[7],p[6],p[5],p[1],p[2],zvec};

		ord[2][1] = OrdV3<F>{p[2],p[6],p[5],p[1],c+p[1],c+p[5],zvec};
		ord[2][2] = OrdV3<F>{p[0],p[4],p[5],p[1],c+p[1],c+p[5],zvec};
		ord[2][4] = OrdV3<F>{p[3],p[0],p[1],p[2],c+p[4],c+p[5],c+p[1]};
		ord[2][5] = OrdV3<F>{p[7],p[4],p[5],p[6],c+p[0],c+p[1],c+p[5]};

		ord[3][1] = OrdV3<F>{p[1],p[5],p[6],p[2],c+p[2],c+p[6],zvec};
		ord[3][3] = OrdV3<F>{p[3],p[7],p[6],p[2],c+p[2],c+p[6],zvec};
		ord[3][4] = OrdV3<F>{p[0],p[3],p[2],p[1],c+p[7],c+p[6],c+p[2]};
		ord[3][5] = OrdV3<F>{p[4],p[7],p[6],p[5],c+p[3],c+p[2],c+p[6]};

		ord[4][1] = OrdV3<F>{p[5],p[6],p[2],p[1],c+p[1],c+p[2],zvec};
		ord[4][2] = OrdV3<F>{p[4],p[5],p[1],p[0],c+p[1],c+p[2],c+p[3]};
		ord[4][3] = OrdV3<F>{p[7],p[6],p[2],p[3],c+p[2],c+p[1],c+p[0]};
		ord[4][4] = OrdV3<F>{p[0],p[3],p[2],p[1],c+p[1],c+p[2],zvec};

		ord[5][1] = OrdV3<F>{p[1],p[2],p[6],p[5],c+p[5],c+p[6],zvec};
		ord[5][2] = OrdV3<F>{p[0],p[4],p[5],p[1],c+p[7],c+p[6],c+p[5]};
		ord[5][3] = OrdV3<F>{p[3],p[7],p[6],p[2],c+p[4],c+p[5],c+p[6]};
		ord[5][5] = OrdV3<F>{p[4],p[7],p[6],p[5],c+p[5],c+p[6],zvec};
		break;
		case 2: // 1,2,1
		ord[0][0] = OrdV3<F>{p[0],p[4],p[7],p[3],c+p[3],c+p[7],zvec};
		ord[0][3] = OrdV3<F>{p[2],p[6],p[7],p[3],c+p[3],c+p[7],zvec};
		ord[0][4] = OrdV3<F>{p[1],p[2],p[3],p[0],c+p[3],c+p[7],c+p[4]};
		ord[0][5] = OrdV3<F>{p[5],p[6],p[7],p[4],c+p[7],c+p[3],c+p[0]};

		ord[1][1] = OrdV3<F>{p[1],p[5],p[6],p[2],c+p[2],c+p[6],zvec};
		ord[1][3] = OrdV3<F>{p[3],p[7],p[6],p[2],c+p[2],c+p[6],zvec};
		ord[1][4] = OrdV3<F>{p[0],p[3],p[2],p[1],c+p[2],c+p[6],c+p[5]};
		ord[1][5] = OrdV3<F>{p[4],p[7],p[6],p[5],c+p[6],c+p[2],c+p[1]};

		ord[2][0] = OrdV3<F>{p[0],p[4],p[7],p[3],p[2],p[6],zvec};
		ord[2][1] = OrdV3<F>{p[1],p[5],p[6],p[2],p[3],p[7],zvec};
		ord[2][3] = OrdV3<F>{p[3],p[2],p[6],p[7],zvec,zvec,zvec};
		ord[2][4] = OrdV3<F>{p[0],p[1],p[2],p[3],p[7],p[6],zvec};
		ord[2][5] = OrdV3<F>{p[4],p[5],p[6],p[7],p[3],p[2],zvec};

		ord[4][0] = OrdV3<F>{p[4],p[7],p[3],p[0],c+p[3],c+p[2],c+p[1]};
		ord[4][1] = OrdV3<F>{p[5],p[6],p[2],p[1],c+p[2],c+p[3],c+p[0]};
		ord[4][3] = OrdV3<F>{p[7],p[6],p[2],p[3],c+p[3],c+p[2],zvec};
		ord[4][4] = OrdV3<F>{p[0],p[1],p[2],p[3],c+p[3],c+p[2],zvec};

		ord[5][0] = OrdV3<F>{p[0],p[3],p[7],p[4],c+p[7],c+p[6],c+p[5]};
		ord[5][1] = OrdV3<F>{p[1],p[2],p[6],p[5],c+p[6],c+p[7],c+p[4]};
		ord[5][3] = OrdV3<F>{p[3],p[2],p[6],p[7],c+p[7],c+p[6],zvec};
		ord[5][5] = OrdV3<F>{p[4],p[5],p[6],p[7],c+p[7],c+p[6],zvec};
		break;
		case 3: // 2,2,1
		ord[0][1] = OrdV3<F>{p[1],p[5],p[6],p[2],c+p[3],c+p[7],zvec};
		ord[0][3] = OrdV3<F>{p[3],p[7],p[6],p[2],c+p[3],c+p[7],zvec};
		ord[0][4] = OrdV3<F>{p[0],p[3],p[2],p[1],c+p[3],c+p[7],c+p[4]};
		ord[0][5] = OrdV3<F>{p[4],p[7],p[6],p[5],c+p[7],c+p[3],c+p[0]};

		ord[2][1] = OrdV3<F>{p[1],p[5],p[6],p[2],c+p[1],c+p[5],zvec};
		ord[2][3] = OrdV3<F>{p[3],p[7],p[6],p[2],c+p[1],c+p[5],zvec};
		ord[2][4] = OrdV3<F>{p[0],p[1],p[2],p[3],c+p[1],c+p[5],c+p[4]};
		ord[2][5] = OrdV3<F>{p[4],p[5],p[6],p[7],c+p[5],c+p[1],c+p[0]};

		ord[4][1] = OrdV3<F>{p[5],p[6],p[2],p[1],c+p[3],c+p[2],c+p[1]};
		ord[4][3] = OrdV3<F>{p[7],p[6],p[2],p[3],c+p[1],c+p[2],c+p[3]};
		ord[4][4] = OrdV3<F>{p[0],p[1],p[2],p[3],c+p[1],c+p[2],c+p[3]};

		ord[5][1] = OrdV3<F>{p[1],p[2],p[6],p[5],c+p[7],c+p[6],c+p[5]};
		ord[5][3] = OrdV3<F>{p[3],p[2],p[6],p[7],c+p[5],c+p[6],c+p[7]};
		ord[5][5] = OrdV3<F>{p[4],p[5],p[6],p[7],c+p[5],c+p[6],c+p[7]};
		break;
		case 4: // 1,1,2
		ord[0][0] = OrdV3<F>{p[0],p[3],p[7],p[4],c+p[4],c+p[7],zvec};
		ord[0][2] = OrdV3<F>{p[1],p[5],p[4],p[0],c+p[4],c+p[7],c+p[3]};
		ord[0][3] = OrdV3<F>{p[2],p[6],p[7],p[3],c+p[7],c+p[4],c+p[0]};
		ord[0][5] = OrdV3<F>{p[5],p[6],p[7],p[4],c+p[4],c+p[7],zvec};

		ord[1][1] = OrdV3<F>{p[1],p[2],p[6],p[5],c+p[5],c+p[6],zvec};
		ord[1][2] = OrdV3<F>{p[0],p[4],p[5],p[1],c+p[5],c+p[6],c+p[2]};
		ord[1][3] = OrdV3<F>{p[3],p[7],p[6],p[2],c+p[6],c+p[5],c+p[1]};
		ord[1][5] = OrdV3<F>{p[4],p[7],p[6],p[5],c+p[5],c+p[6],zvec};

		ord[2][0] = OrdV3<F>{p[3],p[7],p[4],p[0],c+p[4],c+p[5],c+p[1]};
		ord[2][1] = OrdV3<F>{p[2],p[6],p[5],p[1],c+p[5],c+p[4],c+p[0]};
		ord[2][2] = OrdV3<F>{p[0],p[1],p[5],p[4],c+p[4],c+p[5],zvec};
		ord[2][5] = OrdV3<F>{p[7],p[6],p[5],p[4],c+p[4],c+p[5],zvec};

		ord[3][0] = OrdV3<F>{p[0],p[4],p[7],p[3],c+p[7],c+p[6],c+p[2]};
		ord[3][1] = OrdV3<F>{p[1],p[5],p[6],p[2],c+p[6],c+p[7],c+p[3]};
		ord[3][3] = OrdV3<F>{p[3],p[2],p[6],p[7],c+p[7],c+p[6],zvec};
		ord[3][5] = OrdV3<F>{p[4],p[5],p[6],p[7],c+p[7],c+p[6],zvec};

		ord[4][0] = OrdV3<F>{p[0],p[3],p[7],p[4],p[5],p[6],zvec};
		ord[4][1] = OrdV3<F>{p[1],p[2],p[6],p[5],p[4],p[7],zvec};
		ord[4][2] = OrdV3<F>{p[0],p[1],p[5],p[4],p[7],p[6],zvec};
		ord[4][3] = OrdV3<F>{p[3],p[2],p[6],p[7],p[4],p[5],zvec};
		ord[4][5] = OrdV3<F>{p[4],p[5],p[6],p[7],zvec,zvec,zvec};
		break;
		case 5: // 2,1,2
		ord[0][1] = OrdV3<F>{p[1],p[2],p[6],p[5],c+p[4],c+p[7],zvec};
		ord[0][2] = OrdV3<F>{p[0],p[4],p[5],p[1],c+p[4],c+p[7],c+p[3]};
		ord[0][3] = OrdV3<F>{p[3],p[7],p[6],p[2],c+p[7],c+p[4],c+p[0]};
		ord[0][5] = OrdV3<F>{p[4],p[7],p[6],p[5],c+p[4],c+p[7],zvec};

		ord[2][1] = OrdV3<F>{p[2],p[6],p[5],p[1],c+p[4],c+p[5],c+p[1]};
		ord[2][2] = OrdV3<F>{p[0],p[1],p[5],p[4],c+p[1],c+p[5],c+p[4]};
		ord[2][5] = OrdV3<F>{p[7],p[6],p[5],p[4],c+p[1],c+p[5],c+p[4]};

		ord[3][1] = OrdV3<F>{p[1],p[5],p[6],p[2],c+p[7],c+p[6],c+p[2]};
		ord[3][3] = OrdV3<F>{p[3],p[2],p[6],p[7],c+p[2],c+p[6],c+p[7]};
		ord[3][5] = OrdV3<F>{p[4],p[5],p[6],p[7],c+p[2],c+p[6],c+p[7]};

		ord[4][1] = OrdV3<F>{p[1],p[2],p[6],p[5],c+p[1],c+p[2],zvec};
		ord[4][2] = OrdV3<F>{p[0],p[1],p[5],p[4],c+p[1],c+p[2],c+p[3]};
		ord[4][3] = OrdV3<F>{p[3],p[2],p[6],p[7],c+p[2],c+p[1],c+p[0]};
		ord[4][5] = OrdV3<F>{p[4],p[7],p[6],p[5],c+p[1],c+p[2],zvec};
		break;
		case 6: // 1,2,2
		ord[0][0] = OrdV3<F>{p[0],p[3],p[7],p[4],c+p[3],c+p[7],c+p[4]};
		ord[0][3] = OrdV3<F>{p[2],p[6],p[7],p[3],c+p[4],c+p[7],c+p[3]};
		ord[0][5] = OrdV3<F>{p[5],p[6],p[7],p[4],c+p[3],c+p[7],c+p[4]};

		ord[1][1] = OrdV3<F>{p[1],p[2],p[6],p[5],c+p[2],c+p[6],c+p[5]};
		ord[1][3] = OrdV3<F>{p[3],p[7],p[6],p[2],c+p[5],c+p[6],c+p[2]};
		ord[1][5] = OrdV3<F>{p[4],p[7],p[6],p[5],c+p[2],c+p[6],c+p[5]};

		ord[2][0] = OrdV3<F>{p[0],p[4],p[7],p[3],c+p[4],c+p[5],c+p[1]};
		ord[2][1] = OrdV3<F>{p[1],p[5],p[6],p[2],c+p[5],c+p[4],c+p[0]};
		ord[2][3] = OrdV3<F>{p[3],p[2],p[6],p[7],c+p[4],c+p[5],zvec};
		ord[2][5] = OrdV3<F>{p[4],p[5],p[6],p[7],c+p[4],c+p[5],zvec};

		ord[4][0] = OrdV3<F>{p[0],p[3],p[7],p[4],c+p[3],c+p[2],c+p[1]};
		ord[4][1] = OrdV3<F>{p[1],p[2],p[6],p[5],c+p[2],c+p[3],c+p[0]};
		ord[4][3] = OrdV3<F>{p[3],p[2],p[6],p[7],c+p[3],c+p[2],zvec};
		ord[4][5] = OrdV3<F>{p[4],p[5],p[6],p[7],c+p[3],c+p[2],zvec};
		break;
		case 7: // 2,2,2
		ord[0][1] = OrdV3<F>{p[1],p[2],p[6],p[5],c+p[3],c+p[7],c+p[4]};
		ord[0][3] = OrdV3<F>{p[3],p[7],p[6],p[2],c+p[4],c+p[7],c+p[3]};
		ord[0][5] = OrdV3<F>{p[4],p[7],p[6],p[5],c+p[3],c+p[7],c+p[4]};

		ord[2][1] = OrdV3<F>{p[1],p[5],p[6],p[2],c+p[4],c+p[5],c+p[1]};
		ord[2][3] = OrdV3<F>{p[3],p[2],p[6],p[7],c+p[1],c+p[5],c+p[4]};
		ord[2][5] = OrdV3<F>{p[4],p[5],p[6],p[7],c+p[1],c+p[5],c+p[4]};

		ord[4][1] = OrdV3<F>{p[1],p[2],p[6],p[5],c+p[3],c+p[2],c+p[1]};
		ord[4][3] = OrdV3<F>{p[3],p[2],p[6],p[7],c+p[1],c+p[2],c+p[3]};
		ord[4][5] = OrdV3<F>{p[4],p[5],p[6],p[7],c+p[1],c+p[2],c+p[3]};
		break;
	}
}


#endif
