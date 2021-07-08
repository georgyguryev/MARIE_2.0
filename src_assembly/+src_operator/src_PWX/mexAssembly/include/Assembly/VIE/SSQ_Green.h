#ifndef SSQ_GREEN_H
#define SSQ_GREEN_H

#include "EMMC.h"
#include "vec3.h"

template <typename F>
struct PoissonGF {
	F R;
	F G;
	constexpr PoissonGF(vec3<F> const& dr) : R(mag(dr)),G(M_1_4PI/R) { }
};

template <typename F>
struct HelmholtzGF {
	F phi;
	complex<F> G;
	constexpr HelmholtzGF(F const k0,PoissonGF<F> const& p) : phi(k0*p.R),G(complex<F>(cos(phi),-sin(phi))*p.G) { }
};

#endif