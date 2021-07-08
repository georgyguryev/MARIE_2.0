#include <iostream>
#include <stdio.h>
#include <memory>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include "mlaa_boost_types.h"
#include "sie_coup_matrix_element.h"
#include "sie_coup_matrix_tfunc.h"
#include "utils_verbose_support.h"
#include "utils_time_measure.h"
#include "utils_msg_logger.h"
#include "utils_constants.h"
#include "quads_interface_1d.h"
#include "quads_interface_2d.h"
#include "geo_scoil_loader.h"
#include "geo_spatial_domain.h"
#include "geo_rm_file_loader.h"
#include "geo_realistic_model.h"
#include "geo_idx_rm.h"
#include "num_aux_vectors.h"
#include "lr_skelap2d_cross_svd.h"
//#include "lr_skelap2d_cross_mvl.h"
#include "lr_skelap2d_cross_aca.h"


#ifdef _OPENMP
#include <omp.h>
#endif

#include "mex.h"


using std::cout;
using std::endl;
using std::unique_ptr;

using namespace marie_core_0_1;

//////////////////////////////////////////////////////////////////////////////////////////



// template function for  full matrix assembly
template <ViePiecewise vpw_CL, VieGreenFunc  vgf_KL_op>
bool do_full_coupling_assembly(const geoScoilLoader * const p_scoil_loader,
                               const RealisticModel * const p_rmodel,
                               double frequency, mxArray * plhs[]) {

    // total number of independent coordinates [x,y,z] and
    // total number of basis elements per one coordinate
    const size_t q0_max = 3;
    const size_t l_max  = taylor_terms_number<vpw_CL>();

    unique_ptr<CouplingMatrixTargetFunction<vpw_CL, vgf_KL_op>> up_coupling_tfunc
            (new CouplingMatrixTargetFunction<vpw_CL, vgf_KL_op>(string("CouplingMatrixElement")));


    up_coupling_tfunc->set_scoil_loader(p_scoil_loader);
    up_coupling_tfunc->set_realistic_model(p_rmodel);

    size_t n_vie_points = 2;
    size_t m_sie_points = 3;

    up_coupling_tfunc->set_vie_quads_order_1d(n_vie_points);
    up_coupling_tfunc->set_sie_quads_order_2d(m_sie_points);
    up_coupling_tfunc->set_frequency(frequency);

    const int nrows = up_coupling_tfunc->nx_rows();
    const int mcols = up_coupling_tfunc->my_cols();
    const size_t Lx = p_rmodel->cp_sdomain()->Lx();
    const size_t My = p_rmodel->cp_sdomain()->My();
    const size_t Nz = p_rmodel->cp_sdomain()->Nz();


    const uint64_vector & cref_idx_vec = p_rmodel->cp_idx_mask()->cref_spatial_idx_1d_vec();
    const uint64_vector & fort_pos_vec = p_rmodel->cp_idx_mask()->fort_pos_vec();
    const size_t idx_size = cref_idx_vec.size();



    cout << "Coupling dimensions are: "
         << nrows << " x " << mcols << endl;

    // Setup Coupling matrix
    boost_cmplx_2d  Zbc(boost::extents[nrows][mcols]);

#pragma omp parallel for
   for (int vx = 0; vx < nrows; ++vx) {
       for (int elm = 0; elm < mcols; ++elm) {
           Zbc[vx][elm] = up_coupling_tfunc->value(vx, elm);
        }
    }

    plhs[0] = mxCreateDoubleMatrix(nrows, mcols, mxCOMPLEX);

    // setup references to return values
    double* Zr_re = mxGetPr(plhs[0]);
    double* Zr_im = mxGetPi(plhs[0]);

    printf(" Finished coupling matrix assembly!\n");
    printf( "Proceed to copying \n");

    for (size_t q0 = 0; q0 < q0_max; ++q0) {
        for (size_t l = 0; l < l_max; ++l) {
            for (size_t ind = 0; ind < idx_size; ++ind) {

                const int vx =  q0 * l_max * idx_size + l *idx_size +
                        fort_pos_vec[ind];

                const int mat_vx =  ind  +  l * idx_size +
                        q0 * l_max  * idx_size;

                for (int elm = 0; elm < mcols; ++elm) {
                  
                    Zr_re[mat_vx + nrows*elm] = Zbc[vx][elm].real();
                    Zr_im[mat_vx + nrows*elm] = Zbc[vx][elm].imag();

                }
            }
        }
    }

    return true;
}




//////////////////////////////////////////////////////////////////////////////////////////

// template fuction for cross coupling
template <ViePiecewise vpw_CL, VieGreenFunc  vgf_KL_op>
bool do_cross_coupling_assembly(BaseCross2D<dcomplex> * const tcross, const geoScoilLoader * const p_scoil_loader,
                                const RealisticModel * const p_rmodel, double frequency, mxArray * plhs[]) {

    // create target function
    unique_ptr<CouplingMatrixTargetFunction<vpw_CL, vgf_KL_op>>  up_coupling_tfunc
            (new CouplingMatrixTargetFunction<vpw_CL, vgf_KL_op>(string("CouplingMatrixElement")));

    up_coupling_tfunc->set_scoil_loader(p_scoil_loader);
    up_coupling_tfunc->set_realistic_model(p_rmodel);

    size_t n_vie_points = 1;
    size_t m_sie_points = 20;
    size_t nrows, mcols, rank;

    up_coupling_tfunc->set_vie_quads_order_1d(n_vie_points);
    up_coupling_tfunc->set_sie_quads_order_2d(m_sie_points);
    up_coupling_tfunc->set_frequency(frequency);


    const double cross_eps  = 1e-3;

    cout << "Boss, I'm on duty! Sip your coffee and relax!" << endl;

    // run cross coupling
    tcross->set_target_function(up_coupling_tfunc.get());
    tcross->set_truncation_accuracy(cross_eps);
    tcross->do_decompose();

    // get coupling dimensions
    rank  = tcross->rank();
    nrows = tcross->nx_rows();
    mcols = tcross->my_cols();

    cout << "Coupling dimensions are: "
         << nrows << " x " << mcols << endl;

    // allocate Ur and Vr
    plhs[0] = mxCreateDoubleMatrix(nrows, rank, mxCOMPLEX);
    plhs[1] = mxCreateDoubleMatrix(mcols, rank, mxCOMPLEX);

    // setup references to return values
    double* Ur_re = mxGetPr(plhs[0]);
    double* Ur_im = mxGetPi(plhs[0]);
    double* Vr_re = mxGetPr(plhs[1]);
    double* Vr_im = mxGetPi(plhs[1]);

    const dcomplex * const pUr = tcross->ptr_const_U();
    const dcomplex * const pVr = tcross->ptr_const_V();

    // return values
    for (size_t j=0; j < rank; ++j) {
        const double sigma_sqrt =  sqrt(tcross->Sigma(j));
        for (size_t i=0; i < nrows; ++i) {
            Ur_re[i + j * nrows] = sigma_sqrt * pUr[i + j * nrows].real();
            Ur_im[i + j * nrows] = sigma_sqrt * pUr[i + j * nrows].imag();
        }
    }

    for (size_t j=0; j < rank; ++j) {
        const double sigma_sqrt =  sqrt(tcross->Sigma(j));
        for (size_t i=0; i < mcols; ++i) {
            Vr_re[i + j * mcols] = sigma_sqrt * pVr[i + j * mcols].real();
            Vr_im[i + j * mcols] = sigma_sqrt * pVr[i + j * mcols].imag();
        }
    }

    return true;

}

//////////////////////////////////////////////////////////////////////////////////////////

// Gateway function
void  mexFunction(int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]) {

    // set utilities
    Verbose::DebugMemory::setShowFlag(true);
    Verbose::UserSupport::setShowFlag(true);
    Verbose::TechnicalDetails::setShowFlag(true);
    Verbose::TimeStamps::setShowFlag(false);

    TimeMeasure::set_default_stamp_format(TmFormatType::yyyymmdd_hhmmss);
    TimeMeasure::set_default_ssdiff_format(TmFormatType::mmssms);
    TimeMeasure::set_time_separator(":");

    world_timer().push_start_point();
    world_timer().fix_start_point();

    const int nargs = 6;

    if (nargs != nrhs) {
         msg_logger(nothingIdObject(),"Coupling Assembly(): Invalid number of input parameters!");
         return;
    }

    const int number_of_basises   = 2;
    const int number_of_assembly_types = 2;

    // declare
    string body_fname, scoil_fname, loads_fname;
    int  operator_t, basis_type, cross_type;
    double frequency;


    // get model filenames
    body_fname  = string(mxArrayToString(prhs[0])) + ".txt";
    scoil_fname = string(mxArrayToString(prhs[1])) + ".msh";
    loads_fname = string(mxArrayToString(prhs[1])) + ".txt";

    // get basis and cross algorithm settings
    operator_t =  static_cast<int>(mxGetScalar(prhs[2]));
    basis_type =  static_cast<int>(mxGetScalar(prhs[3]));
    cross_type =  static_cast<int>(mxGetScalar(prhs[4]));

    // get frequency
    frequency  = mxGetScalar(prhs[5]);


    cout << "The input  filename is: "<< scoil_fname << endl;
    cout.flush();

    // declare a static Scoil loader
    static unique_ptr <geoScoilLoader> up_scoil_loader;
    static unique_ptr<RealisticModel> up_rmodel;


    // check if the object is initialized first time
    if (nullptr == up_scoil_loader.get()) {
        cout << "Initialize up_scoil_loader" << endl;
        up_scoil_loader.reset(new geoScoilLoader("ScoilLoader"));
        up_scoil_loader->import_ports_properties(loads_fname);
        up_scoil_loader->import_scoil(scoil_fname);
        cout << "Finish initialization \n";
    }

    // check if the object is initialized first time
    if (nullptr == up_rmodel.get()) {
        up_rmodel.reset(new RealisticModel());
        up_rmodel->load_geometry_physical_properties(body_fname, geoRModel_rwMode::FullText);
        up_rmodel->set_frequency(frequency);
        up_rmodel->after_update();
        up_rmodel->cp_idx_mask()->create_fort_pos_for_cidx();
    }

    unique_ptr<BaseCross2D<dcomplex>>    up_scross2d(nullptr);

    // summarize type of assembly
    int assembly_type = basis_type + number_of_basises * cross_type +
                        number_of_basises * number_of_assembly_types * operator_t;

    // run simulations
    switch (assembly_type) {
    case 0:
        cout << "PWC, L_type, Full Assembly" << endl;
        do_full_coupling_assembly<ViePiecewise::Constant, VieGreenFunc::L_type>
                (up_scoil_loader.get(), up_rmodel.get(), frequency, plhs);
        break;
    // PWL & Full
    case 1:
        cout << "PWL, L_type,  Full Assembly" << endl;
        do_full_coupling_assembly<ViePiecewise::Linear, VieGreenFunc::L_type>
                (up_scoil_loader.get(), up_rmodel.get(), frequency, plhs);
        break;
    // PWC & ACA
    case 2:
        cout << "PWC, L_type, ACA Assembly" << endl;
        up_scross2d.reset(new ACA_Cross2D<dcomplex>());
        do_cross_coupling_assembly<ViePiecewise::Constant, VieGreenFunc::L_type>
                (up_scross2d.get(), up_scoil_loader.get(),
                  up_rmodel.get(), frequency, plhs);
        break;

    // PWL & ACA
    case 3:
        cout << "PWL, L_type, ACA Assembly" << endl;
        up_scross2d.reset(new ACA_Cross2D<dcomplex>());
        do_cross_coupling_assembly<ViePiecewise::Linear, VieGreenFunc::L_type>
                (up_scross2d.get(), up_scoil_loader.get(),
                 up_rmodel.get(), frequency, plhs);
        break;

        //Kop PWC & Full
    case 4:
        cout << "PWC, K_type, Full Assembly" << endl;
        do_full_coupling_assembly<ViePiecewise::Constant, VieGreenFunc::K_type>
                (up_scoil_loader.get(), up_rmodel.get(), frequency, plhs);
        break;

        //Kop PWL & Full
    case 5:
        cout << "PWL, K_type, Full Assembly" << endl;
        do_full_coupling_assembly<ViePiecewise::Linear, VieGreenFunc::K_type>
                (up_scoil_loader.get(), up_rmodel.get(),
                 frequency, plhs);
        break;

        //Kop PWC & ACA
    case 6:
        cout << "PWC, K_type, ACA Assembly" << endl;
        up_scross2d.reset(new ACA_Cross2D<dcomplex>());
        do_cross_coupling_assembly<ViePiecewise::Constant, VieGreenFunc::K_type>
                (up_scross2d.get(), up_scoil_loader.get(),
                 up_rmodel.get(), frequency, plhs);
        break;

        //Kop PWL & ACA
    case 7:
        cout << "PWL, K_type, ACA Assembly" << endl;
        up_scross2d.reset(new ACA_Cross2D<dcomplex>());
        do_cross_coupling_assembly<ViePiecewise::Linear, VieGreenFunc::K_type>
                (up_scross2d.get(), up_scoil_loader.get(),
                 up_rmodel.get(), frequency,  plhs);
        break;

    // Full matrix assembly
    default:
        break;
    }

}


//End of file
