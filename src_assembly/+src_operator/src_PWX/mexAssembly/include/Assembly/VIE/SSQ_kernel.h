#ifndef SSQ_KERNEL_H
#define SSQ_KERNEL_H

#include "SSQ_Green.h"

template <typename ...F>
std::common_type_t<F...> const (&cat(F const&... x))[sizeof...(F)] {
	std::common_type_t<F...> const X[sizeof...(F)] = {x...};
	std::common_type_t<F...> const (&R)[sizeof...(F)] = X;
	return R;
}

template <typename F>
struct R4D {
	vec3<F> E;
	vec3<F> G;
	vec3<F> R;
	vec3<F> P;
	constexpr R4D(std::array<vec3<F>,4> const& r)  {
		vec3<F> const m[2] = {r[0]-r[1],r[2]-r[3]};
		vec3<F> const p[2] = {r[0]+r[1],r[2]+r[3]};
		E = -m[0]+m[1];
		G = -p[0]+p[1];
		R = m[0]+m[1];
		P = p[0]+p[1];
	}
	constexpr R4D(vec3<F> const r0,vec3<F> const r1,vec3<F> const r2,vec3<F> const r3) : R4D(std::array<vec3<F>,4>{r0,r1,r2,r3}) { }
	constexpr F jacobian(F const,F const) const;
	constexpr F jacobian(tensor<F,2> const&) const;
	constexpr vec3<F> position(F const&,F const&) const;
	constexpr vec3<F> position(tensor<F,2> const&) const;
};

template <typename F>
constexpr F R4D<F>::jacobian(F const U,F const V) const {
	vec3<F> const Ev = E+V*R;
	vec3<F> const Gv = G+U*R;
	return sqrt(magsq(Ev)*magsq(Gv)-magsq(dot(Ev,Gv)))/16;
}
template <typename F> constexpr F R4D<F>::jacobian(tensor<F,2> const& UV) const {
	return jacobian(UV[0],UV[1]);
}

template <typename F> constexpr vec3<F> R4D<F>::position(F const& U,F const& V) const {
	return ((F)0.25)*(P+U*E+V*G+U*V*R);
}
template <typename F> constexpr vec3<F> R4D<F>::position(tensor<F,2> const& UV) const {
	return position(UV[0],UV[1]);
}

template <typename F>
constexpr F scaled_disp(int l,vec3<F> const& dri,vec3<F> const& r,vec3<F> const& c,vec3<F> const& N) {
	return (l == 0) ? 1 : (r[--l]-c[l])*dri[l]+0.5*N[l];
}

template <typename F>
constexpr vec3<F> triplet_normal(vec3<F> const (&p)[3]) {
	vec3<F> const n = cross(p[1]-p[0],p[2]-p[0]);
	// to-do: follow up with Ioannis
	return n / magsq(n);
}

template <typename F>
inline F unit_normal_dot(int const g,vec3<F> const& x) {
	return (-1+2*(g%2))*x[g/2];
}
#pragma omp declare simd uniform(g,l,x)
template <typename F>
inline F unit_normal_scal(int g,int l,F x) {
	return (-1+2*(g%2))*(l == g/2)*x;
}

template <typename F,int type>
struct VIEKernel {
	F const k0;
	vec3<F> const dri;
	int const f;
	int const g;
	R4D<F>  const r[2];
	vec3<F> const c[2];
	constexpr VIEKernel(F const& k0,vec3<F> const& dri,int f,int g,R4D<F> const (&r)[2],vec3<F> const (&c)[2]) : 
		k0(k0),dri(dri),f(f),g(g),r(r),c(c) { }
	template <bool pwl> void kernel(F const w,tensor<F,2,2> const& Upq,
	F (&Ir)[VIE::Ni<pwl>],F (&Ii)[VIE::Ni<pwl>]) const;
};

template <typename F,int Nq>
struct Gauss2DSS {
	vec3<F> hdr;
	vec3<F> c;
	constexpr Gauss2DSS(vec3<F> const& dr,vec3<F> const& c) : hdr(dr/2),c(c) { }
	#pragma omp declare simd uniform(this,M)
	constexpr vec3<F> x(int const M,int const q1,int const q2) const {
		constexpr int r[3][2] = {
			{1,2},{0,2},{0,1}
		};
		vec3<F> p = c;
		p[r[M][0]] += hdr[r[M][0]]*Gauss<F,Nq>::x(q1);
		p[r[M][1]] += hdr[r[M][1]]*Gauss<F,Nq>::x(q2);
		return p;
	}
};

template <typename F,int Nq,int type>
struct NSKernel {
	F const k0;
	vec3<F> const dr;
	vec3<F> const dri;
	Gauss2DSS<F,Nq> const Q[2];
	int const f;
	int const g;
	constexpr NSKernel(F const k0,vec3<F> const& dr,vec3<F> const& dri,Gauss2DSS<F,Nq> const (&Q)[2],int const f,int const g) : 
		k0(k0),dr(dr),dri(dri),Q(Q),f(f),g(g) { }
	template <bool pwl> void operator()(F (&Ir)[VIE::Ni<pwl>],F (&Ii)[VIE::Ni<pwl>],
		bool_t<pwl>,int q1,int q2,int r1,int r2) const;
};

#pragma omp declare simd uniform(v,c1,c2,l)
template <bool pq,typename F> inline F basis(auto const& v,vec3<F> const& c1,vec3<F> const& c2,vec3<F> const& p,int l) {
	int const h = !pq ? v.f : v.g;
	vec3<F> const& c = !pq ? c1 : c2;
	return (l-- == 0) ? 1 : ((p[l]-c[l])*v.dri[l]+unit_normal_scal(h,l,0.5));
}
#pragma omp declare simd uniform(v,l) 
template <bool pq,typename F,int type> inline F basis(VIEKernel<F,type> const& v,vec3<F> const& p,int l) {
	return basis<pq,F>(v,v.c[0],v.c[1],p,l);
}
#pragma omp declare simd uniform(v,l)
template <bool pq,typename F,int Nq,int type> inline F basis(NSKernel<F,Nq,type> const& v,vec3<F> const& p,int l) {
	return basis<pq,F>(v,v.Q[0].c,v.Q[1].c,p,l);
}

template <typename I=int>
struct UnitN {
	I f;
	constexpr UnitN(I f) : f{f} { }
	constexpr I c() const { return f/2; }
	constexpr I sgn() const { return -1+2*(f%2); }
	constexpr I operator[](I d) const { return I(c() == d); }
	template <typename F>
	constexpr F dot(vec3<F> const& v) const { return sgn()*v[c()]; }
};

#pragma omp declare simd uniform(this)
template <typename F,int type>
template <bool pwl>
inline void VIEKernel<F,type>::kernel(F const w,tensor<F,2,2> const& Upq,
F (&Ir)[VIE::Ni<pwl>],F (&Ii)[VIE::Ni<pwl>]) const {
	vec3<F> const p = r[0].position(Upq[0]);
	vec3<F> const q = r[1].position(Upq[1]);
	F const J = r[0].jacobian(Upq[0])*r[1].jacobian(Upq[1]);
	vec3<F> const dR = p-q;
	UnitN<int> N[2] = {UnitN<int>(f),UnitN<int>(g)};
	PoissonGF<F> const   P = PoissonGF<F>(dR);
	HelmholtzGF<F> const H = HelmholtzGF<F>(k0,P);
	complex<F> const Gf = (H.G*complex<F>(1,H.phi)-P.G)/(H.phi*H.phi);
	
	complex<F> G0;
	switch(type) {
		case 0:
			G0 = H.G;
		break;
		case 1:
		case 2:
			G0 = N[1].dot(dR)*Gf;
		break;
		case 4:
			G0 = N[0].dot(dR)*Gf;
		break;
		case 5:
			G0 = N[0].dot(dR)*(Gf-P.G/2);
		break;
		case 3:
		case 6:
		case 7:
			G0 = (P.G-H.G)/(k0*k0);
		break;
	}
	G0 *= w*J;
	for(int b = 0;b < VIE::Ni<pwl>;b++) {
		F const b0 = basis<0>(*this,p,VIE::Il<pwl>[b]);
		F const b1 = basis<1>(*this,q,VIE::Ip<pwl>[b]);
		switch(type) {
			case 3:
			case 5:
				Ir[b] += real(G0);
				Ii[b] += imag(G0);
			break;
			case 0:
			case 4:
				Ir[b] += b0*b1*real(G0);
				Ii[b] += b0*b1*imag(G0);
			break;
			case 2:
			case 6:
				Ir[b] += b0*real(G0);
				Ii[b] += b0*imag(G0);
			break;
			case 1:
			case 7:
				Ir[b] += b1*real(G0);
				Ii[b] += b1*imag(G0);
			break;
		}
	}
}

template <typename F,int type,SType I> inline VIEKernel<F,type> VIEKernelfromOrd(
F const& k0,vec3<F> const& dri,int const f,int const g,vec3<F> const (&ord)[7]) {
	switch(I) {
		case EAO:
		case EAC: {
			R4D<F> const r4D[2] {
				R4D<F>(std::array<vec3<F>,4>{ord[3],ord[2],ord[5],ord[4]}),
				R4D<F>(std::array<vec3<F>,4>{ord[2],ord[3],ord[0],ord[1]})
			};
			vec3<F> const C[2] = {
				(ord[3]+ord[5])/2,
				(ord[0]+ord[2])/2
			};
			VIEKernel<F,type> const v(k0,dri,f,g,r4D,C);
			return v;
		}
		case VAO:
		case VAC: {
			R4D<F> const r4D[2] = {
				R4D<F>(ord[2],ord[4],ord[5],ord[6]),
				R4D<F>(ord[2],ord[3],ord[0],ord[1])
			};
			vec3<F> const C[2] = {
				(ord[2]+ord[5])/2,
				(ord[0]+ord[2])/2
			};
			VIEKernel<F,type> const v(k0,dri,f,g,r4D,C);
			return v;
		}
		case SIT: {
			R4D<F> const r4D(ord[0],ord[1],ord[2],ord[3]);
			R4D<F> const r4D2[2] = {r4D,r4D};
			vec3<F> const C = (ord[0]+ord[2])/2;
			vec3<F> const C2[2] = {C,C};
			VIEKernel<F,type> const v(k0,dri,f,g,r4D2,C2);
			return v;
		}
	}
}

#pragma omp declare simd uniform(this)
template <typename F,int Nq,int type>
template <bool pwl>
inline void NSKernel<F,Nq,type>::operator()(F (&Ir)[VIE::Ni<pwl>],F (&Ii)[VIE::Ni<pwl>],bool_t<pwl>,int q1,int q2,int r1,int r2) const {
	vec3<F> const rq = Q[0].x(f/2,q1,q2);
	vec3<F> const qq = Q[1].x(g/2,r1,r2);

	vec3<F> const dR = rq-qq;
	PoissonGF<F> const   P = PoissonGF<F>(dR);
	HelmholtzGF<F> const H = HelmholtzGF<F>(k0,P);
	complex<F> const Gp = (H.G*complex<F>(1,H.phi)-P.G)/(H.phi*H.phi);
	complex<F> G0;
	F const w = Gauss<F,Nq>::w(q1,q2,r1,r2);

	switch(type) {
		case 0:
			G0 = H.G;
		break;
		case 1:
		case 2:
			G0 = unit_normal_dot(g,dR)*Gp;
		break;
		case 3:
		case 6:
		case 7:
			G0 = (H.G-P.G)/-(k0*k0);
		break;
		case 4:
			G0 = unit_normal_dot(f,dR)*Gp;
		break;
		case 5:
			G0 = unit_normal_dot(f,dR)*(Gp-P.G/2);
		break;
	}
	G0 *= w;
	
	for(int b = 0;b < VIE::Ni<pwl>;b++) {
		F const dT = basis<0>(*this,rq,VIE::Il<pwl>[b]);
		F const dB = basis<1>(*this,qq,VIE::Ip<pwl>[b]);
		switch(type) {
			case 0:
			case 4:
				Ir[b] += dT*dB*real(G0);
				Ii[b] += dT*dB*imag(G0);
			break;
			case 1:
			case 7:
				Ir[b] += dB*real(G0);
				Ii[b] += dB*imag(G0);
			break;
			case 2:
			case 6:
				Ir[b] += dT*real(G0);
				Ii[b] += dT*imag(G0);
			break;
			case 3:
			case 5:
				Ir[b] += real(G0);
				Ii[b] += imag(G0);
			break;
		}
	}
}

#endif