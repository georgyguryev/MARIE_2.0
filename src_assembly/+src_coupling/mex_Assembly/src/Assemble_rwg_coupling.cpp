#define _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <math.h>
#include "Assemble_coupling.h"
#include "mex.h"

#ifdef _OPENMP
#include <omp.h>
#endif

// Author: Georgy Guryev, Cambridge, MA, 2019    

using namespace std;


///////////////////////////////////////////////////////////////////////////
//! define the scaling coefficients for PWL and PWC basis functions
template <size_t pwx_component> double pwx_scaling(double x, double y, double z);

//! specitialize the template 
template<> double pwx_scaling <0>(double x, double y, double z) {return 1.0;};
template<> double pwx_scaling <1>(double x, double y, double z) {return 0.5 * x;};
template<> double pwx_scaling <2>(double x, double y, double z) {return 0.5 * y;};
template<> double pwx_scaling <3>(double x, double y, double z) {return 0.5 * z;};


///////////////////////////////////////////////////////////////////////////

// define constructor for Constant case
template<>
Assemble_rwg_coupling<Basis_type::Constant>::Assemble_rwg_coupling():
    v_pwx_scaling_{nullptr} 
    {
        v_pwx_scaling_[0] = & pwx_scaling<0>;
    }

// define constructor for Linear case
template<>
Assemble_rwg_coupling<Basis_type::Linear>::Assemble_rwg_coupling():
    v_pwx_scaling_{nullptr, nullptr, nullptr, nullptr} 
    {
        v_pwx_scaling_[0] = &pwx_scaling<0>;
        v_pwx_scaling_[1] = &pwx_scaling<1>;
        v_pwx_scaling_[2] = &pwx_scaling<2>;
        v_pwx_scaling_[3] = &pwx_scaling<3>;
    }        

/////////////////////////////////////////////////////////////////////////// 
    
// define main assembly method    
template <Basis_type CL_basis>
void  Assemble_rwg_coupling<CL_basis>::assemble(mxComplexDouble * Zout, double * Scoord, double * r_rwg,
                                                      double * sie_quads, double * vie_quads,
                                                      double res, const size_t N_l, const size_t N_vox, double ko)
{
    const size_t N3 = 3;
    const size_t N2 = 2;
    const size_t Np_sie = sie_quads[0];  
    const size_t Np_vie = vie_quads[0];

    // get scaling constants 
    double mu_0 = 4 * M_PI * 1e-7;
    double c_0  = 299792458;
    double e_0  = 1.0 / c_0 / c_0 / mu_0; 
    dcomplex ce = dcomplex(0.0,1.0) * ko * c_0 * e_0;
    dcomplex em_scaling = - 1.0 / ce / 4.0 / M_PI;
            
    // allocate memory for quadrature points for positive and negative triangles 
    double * r_src_p   = new double[N3 * Np_sie];
    double * r_src_n   = new double[N3 * Np_sie];
    double * rho_sie_p = new double[N3 * Np_sie];
    double * rho_sie_n = new double[N3 * Np_sie];
    
    
    // get coordinates of vertices 
    double r_p[N3] = {r_rwg[0],          r_rwg[1],          r_rwg[2]         };
    double r_n[N3] = {r_rwg[0 + N3],     r_rwg[1 + N3],     r_rwg[2 + N3]    };
    double r_2[N3] = {r_rwg[0 + 2 * N3], r_rwg[1 + 2 * N3], r_rwg[2 + 2 * N3]};
    double r_3[N3] = {r_rwg[0 + 3 * N3], r_rwg[1 + 3 * N3], r_rwg[2 + 3 * N3]};
    
    // get coordinates of quadrature points on the faces of a given triangle pair
    get_source_coords(r_p, r_2, r_3, sie_quads, r_src_p, rho_sie_p, Np_sie, 1);
    get_source_coords(r_n, r_2, r_3, sie_quads, r_src_n, rho_sie_n, Np_sie, -1);   
    
    // get length of the common edge
    double L_edge = compute_edge_length(r_2, r_3);
    
    // loop over the voxels (this part should run in parallel)
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (size_t i_vox = 0; i_vox < N_vox; ++i_vox) {
        
        // allocate space for observer coordinates
        double r_obs_voxel[N3];

        // get center of a current voxel
        double r_vox_center[N3] = {Scoord[0 + i_vox * N3],
                                   Scoord[1 + i_vox * N3],
                                   Scoord[2 + i_vox * N3]};
        
        // loop other quadrature points along x axis
        for (size_t i_vie = 0; i_vie < Np_vie; ++i_vie) {
            
            // define x-component of the observation point 
            r_obs_voxel[0] = r_vox_center[0] + res / 2.0 * vie_quads[i_vie + Np_vie + 1];
            
            // loop other quadrature points along y axis
            for (size_t j_vie = 0; j_vie < Np_vie; ++j_vie) {
                
                // define y-component of the observation point 
                r_obs_voxel[1] = r_vox_center[1] + res / 2.0 * vie_quads[j_vie + Np_vie + 1];
                
                for (size_t k_vie = 0; k_vie < Np_vie; ++k_vie) {
                    
                    // sdefine y-component of the observation point
                    r_obs_voxel[2] = r_vox_center[2] + res / 2.0 * vie_quads[k_vie + Np_vie + 1];
                
                    // allocate memory for the vector of distances 
                    double R_vec [N3];
                            
                    for (size_t i_tr = 0; i_tr < N2; ++i_tr) {
                        
                        // create a pointer for quadrature points on current triangle
                        double * r_src = nullptr;
                        double * rho   = nullptr;
                        
                        // if current triangle is the first tr. in rwg basis
                        if(0 == i_tr) {
                            r_src = r_src_p;
                            rho   = rho_sie_p;
                        }
                        else {
                            r_src = r_src_n;
                            rho   = rho_sie_n; 
                        }
                            
                        // loop over surface quadrature points
                        for (size_t i_sie = 0; i_sie < Np_sie; ++i_sie) {
                            
                            // compute the radius-vector between the source and observer
                            R_vec[0] = r_obs_voxel[0] - r_src[0 + i_sie * N3];
                            R_vec[1] = r_obs_voxel[1] - r_src[1 + i_sie * N3];
                            R_vec[2] = r_obs_voxel[2] - r_src[2 + i_sie * N3];
                            
                            
                            //compute the distance
                            const double R_dist = vec_norm_l2(R_vec);
                            
                            //get current weights
                            const double weight = vie_quads[i_vie + 1] * vie_quads[j_vie + 1] *
                                                  vie_quads[k_vie + 1] * sie_quads[i_sie + 1] / 8.0;
                            
                            
                            const double R2   = R_dist * R_dist;
                            const double R3   = R_dist * R2;
                            const double koR  = ko * R_dist;
                            const double koR2 = koR * koR;
                            
                            const dcomplex kappa =  exp(-dcomplex(0.0,1.0) * koR) / R3;
                            const dcomplex P = (koR2 - 3.0 * dcomplex(0.0, 1.0) * koR - 3.0) / R2;
                            const dcomplex Q = (koR2 -       dcomplex(0.0, 1.0) * koR - 1.0) / P;
                            
                            const dcomplex mul_ct = kappa * P;
                            
                            const dcomplex  Gxx = R_vec[0] * R_vec[0] - Q;
                            const dcomplex  Gxy = R_vec[0] * R_vec[1];
                            const dcomplex  Gxz = R_vec[0] * R_vec[2];
                            
                            const dcomplex  Gyy = R_vec[1] * R_vec[1] - Q;
                            const dcomplex  Gyz = R_vec[1] * R_vec[2];
                            
                            const dcomplex  Gzz = R_vec[2] * R_vec[2] - Q;
                            
                            dcomplex Kernel_1[3];
                            Kernel_1[0] = Gxx * rho[0 + i_sie * N3] +  Gxy * rho[1 + i_sie * N3] + Gxz * rho[2 + i_sie * N3];
                            Kernel_1[1] = Gxy * rho[0 + i_sie * N3] +  Gyy * rho[1 + i_sie * N3] + Gyz * rho[2 + i_sie * N3];
                            Kernel_1[2] = Gxz * rho[0 + i_sie * N3] +  Gyz * rho[1 + i_sie * N3] + Gzz * rho[2 + i_sie * N3];
                            
                            // loop over basis components
                            for (size_t component = 0; component < N_l; ++component) {
                                
                                
                                const double basis_scaling = v_pwx_scaling_[component](vie_quads[i_vie + Np_vie + 1],
                                                                                       vie_quads[j_vie + Np_vie + 1],
                                                                                       vie_quads[k_vie + Np_vie + 1]);
                                
                                dcomplex update_x = weight * mul_ct * basis_scaling * Kernel_1[0] * L_edge * em_scaling;
                                dcomplex update_y = weight * mul_ct * basis_scaling * Kernel_1[1] * L_edge * em_scaling;
                                dcomplex update_z = weight * mul_ct * basis_scaling * Kernel_1[2] * L_edge * em_scaling;
                                
                                // store x components
                                Zout[i_vox + component * N_vox].real  += real(update_x);
                                Zout[i_vox + component * N_vox].imag  += imag(update_x);

                                
                                // store y components
                                Zout[i_vox + component * N_vox + N_vox * N_l].real += real(update_y);
                                Zout[i_vox + component * N_vox + N_vox * N_l].imag += imag(update_y);

                                
                                // store z components
                                Zout[i_vox + component * N_vox + 2 * N_vox * N_l].real += real(update_z);
                                Zout[i_vox + component * N_vox + 2 * N_vox * N_l].imag += imag(update_z);
                               
                                
                            }
                            
                        }
                    }
                    
                }
            }
        }
        
    }
    
    // de-allocate memory
    delete [] r_src_p;
    delete [] r_src_n;
    delete [] rho_sie_p;
    delete [] rho_sie_n;
    
}

///////////////////////////////////////////////////////////////////////////


void get_source_coords(double * r_v, double * r_2, double * r_3, 
                      double * sie_quads, double * z_sie, double * rho,
                      const size_t Np_sie, int sign) {
    
    const size_t N3 = 3;

    // define index shifts for sie quads
    const size_t sie_s1 = Np_sie + 1;
    const size_t sie_s2 = 2 * Np_sie + 1;
    const size_t sie_s3 = 3 * Np_sie + 1;
        
    for (size_t i_sie = 0; i_sie < Np_sie; ++i_sie) {

        z_sie[0 + i_sie * N3] = r_v[0] * sie_quads[i_sie + sie_s1] + r_2[0] * sie_quads[i_sie + sie_s2] + r_3[0] * sie_quads[i_sie + sie_s3];
        z_sie[1 + i_sie * N3] = r_v[1] * sie_quads[i_sie + sie_s1] + r_2[1] * sie_quads[i_sie + sie_s2] + r_3[1] * sie_quads[i_sie + sie_s3];
        z_sie[2 + i_sie * N3] = r_v[2] * sie_quads[i_sie + sie_s1] + r_2[2] * sie_quads[i_sie + sie_s2] + r_3[2] * sie_quads[i_sie + sie_s3];

        // compute vector rho
        rho[0 + i_sie * N3] = sign * (z_sie[0 + i_sie * N3] - r_v[0]);
        rho[1 + i_sie * N3] = sign * (z_sie[1 + i_sie * N3] - r_v[1]);
        rho[2 + i_sie * N3] = sign * (z_sie[2 + i_sie * N3] - r_v[2]);
        
    }
    
}

///////////////////////////////////////////////////////////////////////////
//! function vec_norm_l2 computes the second norm of the 3D vector r
double vec_norm_l2(double * r) { return sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);}

//! function computes the length of vector, formed by r_1 and r_2
double compute_edge_length(double * r_1, double * r_2) {
    
    // the length of acceptable vector
    const size_t N3 = 3;
    
    // get the vector, formed by the two radius-vectors
    double L[N3] = {r_2[0] - r_1[0], r_2[1] - r_1[1], r_2[2] - r_1[2]};
        
    return vec_norm_l2(L);
    
}

///////////////////////////////////////////////////////////////////////////

template class Assemble_rwg_coupling <Basis_type::Constant>;
template class Assemble_rwg_coupling <Basis_type::Linear  >;

///////////////////////////////////////////////////////////////////////////

Base_Assemble_rwg_coupling * Assemble_rwg_coupling_factory::get_assembler(const size_t & N_l) {
    
        switch (N_l) {
            case 1:
                return new Assemble_rwg_coupling <Basis_type::Constant>;
            case 4:
                return new Assemble_rwg_coupling <Basis_type::Linear  >;
            default:
                return nullptr;
        }  
    }

