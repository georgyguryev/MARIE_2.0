#ifndef SSQ_WEIGHTS_H
#define SSQ_WEIGHTS_H

#include "VIE.h"

template <typename F,bool pwl>
struct SS_WK {
	constexpr static int p[3] = {2,0,1};
	constexpr static int q[3] = {1,2,0};
	constexpr static vec3<int> N[6] {
		{-1, 0, 0},{+1, 0, 0},
		{ 0,-1, 0},{ 0,+1, 0},
		{ 0, 0,-1},{ 0, 0,+1}
	};
	constexpr static vec3<int> I[3] {
		{+1,0,0},{0,+1,0},{0,0,+1}
	};
	vec3<F> H[4];
	vec3<int> Nf;
	vec3<int> Ng;
	vec3<int> pp;
	vec3<int> qq;
	F k0;
	F W0;
	F W1[VIE::Ni<pwl>];
	F W2[VIE::Nb<pwl>];
	F W3[VIE::Nb<pwl>];
	bool Z[4];
	SS_WK(vec3<F> const& dr,F const k0,int const f,int const g,int const i) :
		H({{0,0,0},{1/dr.x,0,0},{0,1/dr.y,0},{0,0,1/dr.z}}),
		Nf(N[f]),Ng(N[g]),pp(I[p[i]]),qq(I[q[i]]),k0(k0),W0{},W1{},W2{},W3{} {
		W0 = dot(cross(pp,qq),Ng);
		for(int b = 0;b < VIE::Ni<pwl>;b++) {
			int const l = VIE::Il<pwl>[b];
			int const p = VIE::Ip<pwl>[b];
			W1[b] = dot(cross(qq,pp),H[p])*dot(Ng,H[l]) / -(k0*k0);
		}
		for(int m = 0;m < VIE::Nb<pwl>;m++) {
			W2[m] = dot(Nf,Ng)*dot(H[m],cross(pp,qq));
		}
		vec3<F> const A(
			dot(cross(qq,I[0]),Ng),
			dot(cross(qq,I[1]),Ng),
			dot(cross(qq,I[2]),Ng)
		);
		for(int l = 0;l < VIE::Nb<pwl>;l++) {
			vec3<F> const B = dot(H[l],Nf)*pp;
			W3[l] = dot(A,B);
		}
		Z[0] = W0 == 0;
		Z[1] = true;
		Z[2] = true;
		Z[3] = true;
		for(int b = 0;b < VIE::Ni<pwl> && Z[1];b++) {
			Z[1] = Z[1] && W1[b] == 0;
		}
		for(int m = 0;m < VIE::Nb<pwl> && Z[2];m++) {
			Z[2] = Z[2] && W2[m] == 0;
		}
		for(int l = 0;l < VIE::Nb<pwl> && Z[3];l++) {
			Z[3] = Z[3] && W3[l] == 0;
		}
	}
	template <int type>
	constexpr F operator()(int_t<type>,int b) const {
		switch(type) {
			case 4:
				return W0;
			case 5:
				return W1[b];
			case 6:
				return W2[VIE::Ip<pwl>[b]];
			case 7:
				return W3[VIE::Il<pwl>[b]];
			default:
				return 0;
		}
	}
};
template <typename F,bool pwl> constexpr int SS_WK<F,pwl>::p[3];
template <typename F,bool pwl> constexpr int SS_WK<F,pwl>::q[3];
template <typename F,bool pwl> constexpr vec3<int> SS_WK<F,pwl>::N[6];
template <typename F,bool pwl> constexpr vec3<int> SS_WK<F,pwl>::I[3];

template <typename F,bool pwl>
struct SS_WN {
	constexpr static int p[6] = {0,0,0,1,1,2};
	constexpr static int q[6] = {0,1,2,1,2,2};
	constexpr static vec3<int> N[6] {
		{-1, 0, 0},{+1, 0, 0},
		{ 0,-1, 0},{ 0,+1, 0},
		{ 0, 0,-1},{ 0, 0,+1}
	};
	constexpr static vec3<int> I[3] {
		{+1,0,0},{0,+1,0},{0,0,+1}
	};
	vec3<F> H[4];
	vec3<int> Nf;
	vec3<int> Ng;
	vec3<int> pp;
	vec3<int> qq;
	F W0;
	F W1[VIE::Nb<pwl>];
	F W2[VIE::Nb<pwl>];
	F W3[VIE::Ni<pwl>];
	bool Z[4];
	constexpr SS_WN(vec3<F> const& dr,F const k0,int const f,int const g,int const i) :
		H({{0,0,0},{1/dr.x,0,0},{0,1/dr.y,0},{0,0,1/dr.z}}),
		Nf(N[f]),Ng(N[g]),pp(I[p[i]]),qq(I[q[i]]),W0{},W1{},W2{},W3{} {
		W0 = dot(cross(Nf,pp),cross(Ng,qq));
		for(int l = 0;l < VIE::Nb<pwl>;l++) {
			W1[l] = dot(cross(cross(H[l],pp),qq),Nf);
		}
		for(int m = 0;m < VIE::Nb<pwl>;m++) {
			W2[m] = dot(cross(pp,cross(H[m],qq)),Nf);
		}
		for(int b = 0;b < VIE::Ni<pwl>;b++) {
			int const l = VIE::Il<pwl>[b];
			int const p = VIE::Ip<pwl>[b];
			vec3<F> const a = cross(H[l],pp);
			vec3<F> const Q[3] = {cross(a,I[0]),cross(a,I[1]),cross(a,I[2])};
			vec3<F> const W[3] = {qq.x*H[p],qq.y*H[p],qq.z*H[p]};
			W3[b] = dot(Q[0],Nf)*dot(W[0],Ng)+dot(Q[1],Nf)*dot(W[1],Ng)+dot(Q[2],Nf)*dot(W[2],Ng);
		}
		Z[0] = W0 == 0;
		Z[1] = true;
		Z[2] = true;
		Z[3] = true;
		for(int l = 0;l < VIE::Nb<pwl> && Z[1];l++) {
			Z[1] = Z[1] && W1[l] == 0;
		}
		for(int l = 0;l < VIE::Nb<pwl> && Z[2];l++) {
			Z[2] = Z[2] && W2[l] == 0;
		}
		for(int b = 0;b < VIE::Ni<pwl> && Z[3];b++) {
			Z[3] = Z[3] && W3[b] == 0;
		}
	}
	template <int type>
	constexpr F operator()(int_t<type>,int b) const {
		switch(type) {
			case 0:
				return W0;
			case 1:
				return W1[VIE::Il<pwl>[b]];
			case 2:
				return W2[VIE::Ip<pwl>[b]];
			case 3:
				return W3[b];
			default:
				return 0;
		}
	}
};
template <typename F,bool pwl> constexpr int SS_WN<F,pwl>::p[6];
template <typename F,bool pwl> constexpr int SS_WN<F,pwl>::q[6];
template <typename F,bool pwl> constexpr vec3<int> SS_WN<F,pwl>::N[6];
template <typename F,bool pwl> constexpr vec3<int> SS_WN<F,pwl>::I[3];

template <typename F,bool pwl,VIE::OpType O>
using SS_WO = std::conditional_t<VIE::isN<O>,SS_WN<F,pwl>,SS_WK<F,pwl>>;
template <typename F,bool pwl,int KType>
using SS_WI = SS_WO<F,pwl,KType < 4 ? VIE::N : VIE::K>;

#endif