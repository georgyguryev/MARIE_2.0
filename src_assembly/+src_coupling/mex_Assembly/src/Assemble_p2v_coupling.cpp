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
Assemble_p2v_coupling<Basis_type::Constant>::Assemble_p2v_coupling():
    v_pwx_scaling_{nullptr} 
    {
        v_pwx_scaling_[0] = & pwx_scaling<0>;
    }

// define constructor for Linear case
template<>
Assemble_p2v_coupling<Basis_type::Linear>::Assemble_p2v_coupling():
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
void  Assemble_p2v_coupling<CL_basis>::assemble(mxComplexDouble * Z_Nop_out, mxComplexDouble * Z_Kop_out, double * Scoord,
                                                double * r_points, double * vie_quads,
                                                double res, const size_t N_l, const size_t N_vox,
                                                const size_t N_pts, double ko)
{
    
    const size_t N3 = 3;
    const size_t Np_vie  = vie_quads[0];
    const size_t N_vox3  = N3 * N_vox;
    const size_t N_vox_l = N_vox * N_l;
    const size_t N_dofs  = N_l * N_vox3; 
    
    
    // get scaling constants 
    double mu_0 = 4 * M_PI * 1e-7;
    double c_0  = 299792458;
    double e_0  = 1.0 / c_0 / c_0 / mu_0; 
    dcomplex ce = dcomplex(0.0,1.0) * ko * c_0 * e_0;
    dcomplex e_scaling = - 1.0 / ce / 4.0 / M_PI;
    dcomplex m_scaling = - 1.0 / 4.0 / M_PI;

    for (size_t i_point = 0; i_point < N_pts; ++i_point) {
        
        double r_src[N3] = {r_points[i_point], 
                            r_points[i_point + N_pts],
                            r_points[i_point + 2* N_pts]};
    
        // loop over the voxels (this part should run in parallel)
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
                       
                        
                        // compute the radius-vector between the source and observer
                        R_vec[0] = r_obs_voxel[0] - r_src[0];
                        R_vec[1] = r_obs_voxel[1] - r_src[1];
                        R_vec[2] = r_obs_voxel[2] - r_src[2];
                        
                        
                        //compute the distance
                        const double R_dist = vec_norm_l2(R_vec);
                        
                        //get current weights
                        const double weight = vie_quads[i_vie + 1] * vie_quads[j_vie + 1] *
                                vie_quads[k_vie + 1] / 8.0;
                        
                        
                        const double R2   = R_dist * R_dist;
                        const double R3   = R_dist * R2;
                        const double koR  = ko * R_dist;
                        const double koR2 = koR * koR;
                        
                        const dcomplex kappa =  exp(-dcomplex(0.0,1.0) * koR) / R3;
                        const dcomplex gamma = dcomplex(0.0, 1.0) * koR + 1.0;
                        const dcomplex P = (koR2 - 3.0 * dcomplex(0.0, 1.0) * koR - 3.0) / R2;
                        const dcomplex Q = (koR2 -       dcomplex(0.0, 1.0) * koR - 1.0) / P;
                        
                        const dcomplex mul_ct = kappa * P;
                        
                        // define dyadic Green's function for K operator
                        const dcomplex G = kappa * gamma * m_scaling;
                        
                        const dcomplex  Gxx = R_vec[0] * R_vec[0] - Q;
                        const dcomplex  Gxy = R_vec[0] * R_vec[1];
                        const dcomplex  Gxz = R_vec[0] * R_vec[2];
                        
                        const dcomplex  Gyy = R_vec[1] * R_vec[1] - Q;
                        const dcomplex  Gyz = R_vec[1] * R_vec[2];
                        
                        const dcomplex  Gzz = R_vec[2] * R_vec[2] - Q;
                        
                        // loop over basis components
                        for (size_t component = 0; component < N_l; ++component) {
                            
                            
                            const double basis_scaling = v_pwx_scaling_[component](vie_quads[i_vie + Np_vie + 1],
                                                                                   vie_quads[j_vie + Np_vie + 1],
                                                                                   vie_quads[k_vie + Np_vie + 1]);
                            
                            
                            dcomplex E_Nop_xx = weight * mul_ct * basis_scaling * Gxx * e_scaling;
                            dcomplex E_Nop_xy = weight * mul_ct * basis_scaling * Gxy * e_scaling;
                            dcomplex E_Nop_xz = weight * mul_ct * basis_scaling * Gxz * e_scaling;
                            dcomplex H_Kop_yx = weight * basis_scaling * G * R_vec[2];
                            dcomplex H_Kop_xz = weight * basis_scaling * G * R_vec[1];
                            dcomplex H_Kop_zy = weight * basis_scaling * G * R_vec[0];
                            
                                                        
                            dcomplex E_Nop_yy = weight * mul_ct * basis_scaling * Gyy * e_scaling;
                            dcomplex E_Nop_yz = weight * mul_ct * basis_scaling * Gyz * e_scaling;
                            dcomplex E_Nop_zz = weight * mul_ct * basis_scaling * Gzz * e_scaling;
                            
                            size_t idx_x = i_vox  + component * N_vox + i_point * N_dofs * N3;
                            size_t idx_y = i_vox  + component * N_vox + (i_point * N3 + 1) * N_dofs;
                            size_t idx_z = i_vox  + component * N_vox + (i_point * N3 + 2) * N_dofs;
                            
                            
                            // store x components:
                            
                            // E_x and H_x due to x
                            Z_Nop_out[idx_x].real  += real(E_Nop_xx);
                            Z_Nop_out[idx_x].imag  += imag(E_Nop_xx);
                            Z_Kop_out[idx_x].real  += 0.0;
                            Z_Kop_out[idx_x].imag  += 0.0;
                            
                            // E_x due to y
                            Z_Nop_out[idx_y].real  += real(E_Nop_xy);
                            Z_Nop_out[idx_y].imag  += imag(E_Nop_xy);
                            Z_Kop_out[idx_y].real  += -real(H_Kop_yx);
                            Z_Kop_out[idx_y].imag  += -imag(H_Kop_yx);
                            
                            
                            // E_x due to z
                            Z_Nop_out[idx_z].real  += real(E_Nop_xz);
                            Z_Nop_out[idx_z].imag  += imag(E_Nop_xz);
                            Z_Kop_out[idx_z].real  += real(H_Kop_xz);
                            Z_Kop_out[idx_z].imag  += imag(H_Kop_xz);
                            
                            // store y components:
                            
                            //  E_y due to x
                            Z_Nop_out[idx_x + N_vox_l].real  += real(E_Nop_xy);
                            Z_Nop_out[idx_x + N_vox_l].imag  += imag(E_Nop_xy);
                            Z_Kop_out[idx_x + N_vox_l].real  += real(H_Kop_yx);
                            Z_Kop_out[idx_x + N_vox_l].imag  += imag(H_Kop_yx);
                            
                            // E_y due to y
                            Z_Nop_out[idx_y + N_vox_l].real  += real(E_Nop_yy);
                            Z_Nop_out[idx_y + N_vox_l].imag  += imag(E_Nop_yy);
                            Z_Kop_out[idx_y + N_vox_l].real  += 0.0;
                            Z_Kop_out[idx_y + N_vox_l].imag  += 0.0;
                            
                            // E_y due to z
                            Z_Nop_out[idx_z + N_vox_l].real  += real(E_Nop_yz);
                            Z_Nop_out[idx_z + N_vox_l].imag  += imag(E_Nop_yz);
                            Z_Kop_out[idx_z + N_vox_l].real  += -real(H_Kop_zy);
                            Z_Kop_out[idx_z + N_vox_l].imag  += -imag(H_Kop_zy);
                            
                            // store z components:
                           
                            // E_z due to x
                            Z_Nop_out[idx_x + 2 * N_vox_l].real  += real(E_Nop_xz);
                            Z_Nop_out[idx_x + 2 * N_vox_l].imag  += imag(E_Nop_xz);
                            Z_Kop_out[idx_x + 2 * N_vox_l].real  += -real(H_Kop_xz);
                            Z_Kop_out[idx_x + 2 * N_vox_l].imag  += -imag(H_Kop_xz);
                            
                            
                            // E_z due to y
                            Z_Nop_out[idx_y + 2 * N_vox_l].real  += real(E_Nop_yz);
                            Z_Nop_out[idx_y + 2 * N_vox_l].imag  += imag(E_Nop_yz);
                            Z_Kop_out[idx_y + 2 * N_vox_l].real  += real(H_Kop_zy);
                            Z_Kop_out[idx_y + 2 * N_vox_l].imag  += imag(H_Kop_zy);
                            
                            
                            // E_z due to z
                            Z_Nop_out[idx_z + 2 * N_vox_l].real  += real(E_Nop_zz);
                            Z_Nop_out[idx_z + 2 * N_vox_l].imag  += imag(E_Nop_zz);
                            Z_Kop_out[idx_z + 2 * N_vox_l].real  += 0.0;
                            Z_Kop_out[idx_z + 2 * N_vox_l].imag  += 0.0;
                            
                        }
                        
                    }
                }
                
            }
        }
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

template class Assemble_p2v_coupling <Basis_type::Constant>;
template class Assemble_p2v_coupling <Basis_type::Linear  >;

///////////////////////////////////////////////////////////////////////////

Base_Assemble_p2v_coupling * Assemble_p2v_coupling_factory::get_assembler(const size_t & N_l) {
    
        switch (N_l) {
            case 1:
                return new Assemble_p2v_coupling <Basis_type::Constant>;
            case 4:
                return new Assemble_p2v_coupling <Basis_type::Linear  >;
            default:
                return nullptr;
        }  
    }