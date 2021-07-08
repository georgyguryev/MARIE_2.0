#include <iostream>
#include <complex>
#include "mex.h"

// Author: Georgy Guryev, Cambridge, MA, 2019    

// aliases for data types
using dcomplex = std::complex <double>;
using fcomplex = std::complex <float>;


// define enum for 
enum class Basis_type: int{Constant, Linear, undefined};

// define a function pointer 
using ptr_pwx_components = double (*) (double , double , double );

///////////////////////////////////////////////////////////////////////////

//! define a template that returns the number of dofs per direction
template< Basis_type CL_basis> constexpr size_t  get_num_of_basis_components() noexcept;

// specialization of a template
template<> constexpr size_t  get_num_of_basis_components <Basis_type::Constant>() noexcept {return 1;}
template<> constexpr size_t  get_num_of_basis_components <Basis_type::Linear  >() noexcept {return 4;}



///////////////////////////////////////////////////////////////////////////

//!Base class for RWG 2 voxel coupling

class Base_Assemble_rwg_coupling {
    
public:
    
    // constructor of virtual base class
    Base_Assemble_rwg_coupling(){}
    
    // signature of virtual/interface method
    virtual void assemble(mxComplexDouble * Zout, double * Scoord, double * r_rwg,
                         double * sie_quads, double * vie_quads,
                         double res, const size_t N_l, const size_t N_vox, double ko){}
    
};

//-------------------------------------------------------------------------------

class Base_Assemble_rwg_coupling_matrix {
    
public:
    
    // constructor of virtual base class
    Base_Assemble_rwg_coupling_matrix(){}
    
    // signature of virtual/interface method
    virtual void assemble(mxComplexDouble * Zout, mxComplexDouble * Zbc_Kop,
                          double * Scoord, double * r_p, double * r_n, double * r_2, 
                          double * r_3, double * sie_quads, double * vie_quads, 
                          double res, const size_t N_l, const size_t N_vox,
                          const size_t N_rwg, double ko){}
    
};

//-------------------------------------------------------------------------------

//! Base class for point 2 voxel coupling (used for computation of basis)

class Base_Assemble_p2v_coupling {

public:
    
    // empty constructor 
    Base_Assemble_p2v_coupling(){};
    
    // signature of virtual method for coupling assembly
    virtual void assemble(mxComplexDouble * Z_Nop_out, mxComplexDouble * Z_Kop_out, double * Scoord,
                          double * r_points, double * vie_quads,
                          double res, const size_t N_l, const size_t N_vox, 
                          const size_t N_pts, double ko) {}
   
};


//----------------------------------------------------------------------//

//! Template class for RWG 2 voxel coupling
template <Basis_type CL_basis> 
class Assemble_rwg_coupling: public Base_Assemble_rwg_coupling {
    
public:
    
    // constructor of a derived template class
    Assemble_rwg_coupling();
    
    // main function for rwg coupling assembly
    void  assemble(mxComplexDouble * Zout, double * Scoord, double * r_rwg,
                         double * sie_quads, double * vie_quads,
                         double res, const size_t N_l, const size_t N_vox, double ko);

private:
    
    // define a vector of pointers to the template functions
    ptr_pwx_components v_pwx_scaling_[get_num_of_basis_components<CL_basis>()];
                    
};

//----------------------------------------------------------------------//

//! Template class for RWG 2 voxel coupling
template <Basis_type CL_basis> 
class Assemble_rwg_coupling_matrix: public Base_Assemble_rwg_coupling_matrix {
    
public:
    
    // constructor of a derived template class
    Assemble_rwg_coupling_matrix();
    
    // main function for rwg coupling assembly
    void  assemble(mxComplexDouble * Zout, mxComplexDouble * Zbc_Kop,
                          double * Scoord, double * r_p, double * r_n, double * r_2, 
                          double * r_3, double * sie_quads, double * vie_quads, 
                          double res, const size_t N_l, const size_t N_vox,
                          const size_t N_rwg, double ko);

private:
    
    // define a vector of pointers to the template functions
    ptr_pwx_components v_pwx_scaling_[get_num_of_basis_components<CL_basis>()];
                    
};
//--------------------------------------------------------------------------

//! Template class for point 2 voxel coupling
template<Basis_type CL_basis>
class Assemble_p2v_coupling: public Base_Assemble_p2v_coupling {
    
public:
    
    // constructor to the template class 
    Assemble_p2v_coupling();
    
    //main function for point 2 voxel coupling assembly
    virtual void assemble(mxComplexDouble * Z_Nop_out, mxComplexDouble * Z_Kop_out, double * Scoord,
                          double * r_points, double * vie_quads,
                          double res, const size_t N_l, const size_t N_vox, 
                          const size_t N_pts, double ko);
                          
private:

    // define a vector of pointers to the template functions
    ptr_pwx_components v_pwx_scaling_[get_num_of_basis_components<CL_basis>()];
    
};


///////////////////////////////////////////////////////////////////////////

// declaration of factory class (method get_assembler() creates an appropriate type of coupling)
class Assemble_rwg_coupling_factory {
    
public:
    
    // constructor of factory class
    Assemble_rwg_coupling_factory(){}
    
    // factory method
    Base_Assemble_rwg_coupling * get_assembler(const size_t & N_l);
    
    
};

///////////////////////////////////////////////////////////////////////////

// declaration of factory class (method get_assembler() creates an appropriate type of coupling)
class Assemble_rwg_coupling_matrix_factory {
    
public:
    
    // constructor of factory class
    Assemble_rwg_coupling_matrix_factory(){}
    
    // factory method
    Base_Assemble_rwg_coupling_matrix * get_assembler(const size_t & N_l);
    
    
};

///////////////////////////////////////////////////////////////////////////

// declaration of factory class (method get_assembler() creates an appropriate type of coupling)
class Assemble_p2v_coupling_factory {
    
public:
    
    // constructor of factory class
    Assemble_p2v_coupling_factory(){}
    
    // factory method
    Base_Assemble_p2v_coupling * get_assembler(const size_t & N_l);
    
    
};
///////////////////////////////////////////////////////////////////////////

//! signature of the function that computes sie quadrature points on the triangle
void get_source_coords(double * r_v, double * r_2, double * r_3, 
                      double * sie_quads, double * z_sie, double * rho,
                      const size_t Np_sie, int sign); 


void get_source_coords_mat(double * r_v, double * r_2, double * r_3, 
                      double * sie_quads, double * z_sie, double * rho,
                      const size_t Np_sie, const size_t N_rwg, int sign);
///////////////////////////////////////////////////////////////////////////

//! auxiliary function to compute second norm of tree-dimensional vector
double vec_norm_l2(double * r);

///////////////////////////////////////////////////////////////////////////

//! function computes the length of vector, formed by r_1 and r_2
double compute_edge_length(double * r_1, double * r_2);

double * compute_edge_length_v(double * r_1, double * r_2, size_t N_rwg);


