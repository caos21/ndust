/*
 * Copyright 2017 <Benjamin Santos> <caos21@gmail.com>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef NEVO_H
#define NEVO_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <functional>
#include <cmath>
#include <tuple>
#include <array>
#include <map>

// hdf5 c++ bindings
#include <H5Cpp.h>

// BOOST odeint
#include <boost/numeric/odeint.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

// BOOST Math toolkit
#include <boost/math/tools/roots.hpp>

// #include <cvode/cvode.h>
// // error, check version
// //#include <cvode/cvode_dense.h>
// #include <cvode/cvode_direct.h>
// #include <cvode/cvode_diag.h>
// #include <nvector/nvector_serial.h>
// #include <sundials/sundials_types.h>
// #include <sundials/sundials_math.h>

#include "../include/log.h"
#include "../include/h5plasma.h"
#include "../include/eta.h"
#include "../include/CRate.h"
#include "../include/GridModel.h"
#include "../include/PlasmaModel.h"
#include "../include/NanoModel.h"

#include "../include/Plasma.h"

typedef double value_type;
typedef std::vector<double> state_type;
namespace odeint = boost::numeric::odeint;
typedef odeint::runge_kutta_cash_karp54< state_type > error_stepper_type;
typedef odeint::controlled_runge_kutta< odeint::runge_kutta_cash_karp54< state_type > > stepper_type;

extern "C" {
#include "common.h"
#include "lsoda.h"
#include "lsoda_internal.h"
#include "blas.h" 
}


// forward declaration of qsystem
struct qsystem;

struct ls_qsystem;

// class Solver;

/**
  * class NEvo
  *
  * Represents the data for nanoparticle growth evolution
  *
  */
class NEvo {
public:
// Constructors/Destructors
//
  //NEvo(): moments_file(NULL), qsys(NULL), ls_qsys(NULL), sol(NULL), sih4_file(NULL) {}
  NEvo(): moments_file(NULL), qsys(NULL), ls_qsys(NULL), sih4_file(NULL) {}
  //! Constructor for NEvo
  /*!
    @param  dirname_ Directory for output files.
    @param  grid_filename_ Prefix for output files.
    @param  lg_ Logger instance.
  */
  // NEvo(std::string dirname_,
  //      std::string grid_filename_,
  //      std::string plasma_filename_,
  //      std::string nano_filename_,
  //       src::severity_logger< severity_level > lg_)
  //       : dirname(dirname_),
  // 	  grid_filename(grid_filename_),
  // 	  plasma_filename(plasma_filename_),
  // 	  nano_filename(nano_filename_),
  // 	  lg(lg_), moments_file(NULL), qsys(NULL),
  // 	  ls_qsys(NULL), sol(NULL), sih4_file(NULL) {
  NEvo(std::string dirname_,
       std::string grid_filename_,
       std::string plasma_filename_,
       std::string nano_filename_,
        src::severity_logger< severity_level > lg_)
        : dirname(dirname_),
	  grid_filename(grid_filename_),
	  plasma_filename(plasma_filename_),
	  nano_filename(nano_filename_),
	  lg(lg_), moments_file(NULL), qsys(NULL),
	  ls_qsys(NULL), sih4_file(NULL) {
//     BOOST_LOG_SEV(lg, debug) << "NEvo instantiation";
  }

// Public methods
//
  //! Open h5 datafile
  /*! 
   * Open the file grid_filename.h5 in dirname
  */
  int open();

  //! Close h5 datafile
  /*! 
   * Close the file grid_filename.h5 in dirname
  */
  int close();

  //! Read h5 datafile
  /*! 
   * Read the file grid_filename.h5 in dirname
  */
  int read();

  //! Start calculations
  /*! 
   * Nanoparticle growth evolution
  */
  int evolve();

  //! Start calculations
  /*! 
   * Nanoparticle growth evolution openmp ready
  */
  int evolve_omp();

  //! Start calculations
  /*! 
   * Nanoparticle growth evolution openmp ready
   * with plasma TESTING
  */
  int evolve_selfconsistent();

  //! Start calculations
  /*! 
   * Nanoparticle growth evolution openmp ready
   * with plasma TESTING
  */
  int evolve_onestepselfconsistent(double ctime);

  //! Start calculations
  /*! 
   * Nanoparticle growth evolution openmp ready
  */
  int evolve_radapt();

  //! Plasma init
  /*! 
   * Plasma init reactions
  */
  int init_plasma();

  //! Plasma test
  /*! 
   * Plasma test
  */
  int test_plasma();

  //! Plasma evolution only
  /*! 
   * Plasma evolution
  */
  int evolve_plasma(double ctime);
  

  //! write h5 datafile
  /*! 
   * Write on file grid_filename.h5 in dirname
  */
  int write();
//
// Public attributes
//
  std::string dirname;                //!< string for directory for output.
  std::string grid_filename;        //!< string for prefix for output.
  std::string plasma_filename;        //!< string for prefix for output.
  std::string nano_filename;        //!< string for prefix for output.
  src::severity_logger< severity_level > lg; //!< Logger instance.
  H5::H5File h5obj;                   //!< HDF5 writable object.

  H5::H5File h5obj_plasma;                   //!< HDF5 writable object.

  H5::H5File h5obj_nano;                   //!< HDF5 writable object.

//   GridModel gm;                       //!< Data model for a grid.
  CRate cr;                           //!< Coagulation rate object.

  PlasmaModel pm;                     //!< Plasma model for calculations.

  NanoModel nm;                       //!< Nano model for calculations.

  Plasma plasma;                      //!< Plasma module.

  double ctime;                       //!< Current time of simulation.

  boost_array2d idens;

  // auxiliar densities
  boost_array2d pdens;
  boost_array2d ndens;
  boost_array2d pdens_aux;
  boost_array2d ndens_aux;

  // total rate
  boost_array2d trate2d;
  boost_array2d dt2d;
  
  // rate of coagulation surface growth charging nucleation
  boost_array2d crate2d;
  boost_array2d srate2d;
  boost_array2d qrate2d;
  boost_array2d nrate2d;

  boost_array2d srate2d_aux;
  boost_array2d nrate2d_aux;
  boost_array2d zero2d;

  // Collision frequencies
  boost_array2d efreq;
  boost_array2d ifreq;

  // tunnel current
  boost_array2d tfreq;

  // nanoparticle potential
  boost_array2d phid;

  double efreqfactor;
  double ifreqfactor;

  darray adim_srate;// surface rate 1/s

  darray surface_rate;// surface rate growth in m3/s
  darray surface_area;// area of nanoparticles

  boost_array2d gsurfacegrowth;
  boost_array2d kcoagulation;
  boost_array2d cfrequency;
  boost_array2d jnucleation;

  // vector of birth of particles in section
  boost_array2d birth_vector;
  // vector of death of particles in section
  boost_array2d death_vector;

  // stream for plain dat moments file
  std::fstream* moments_file;
  
  // stream for plain dat plasma file
  std::fstream* plasma_file;

  darray moments;

//   state_type nqdens;
//   friend struct qsystem;
  qsystem* qsys;

  ls_qsystem* ls_qsys;
  
  //friend class Solver;

  //Solver* sol;

  int check;

  // original nano delta time
  double orig_dtn;  

  // ==================== SiH4 ====================
  double nsih4_ini;
  double nsih4;// SiH4 density
  double nsih4_aux;
  double sih4rate;
  double vsih4;// SiH4 thermal velocity
  double sih4_vol;// SiH4 volume
  double sgrowth_rate_sih4;
  double sgrowth_effective; // effective surface growth coefficient
  double sgrowth_total_rate;
  
  std::fstream* sih4_file;
//
private:
// Private methods
//

  //! One step evolution
  /*! 
   * Nanoparticle growth one step
  */
  int evolve_one_step(double ctime);

  //! One step evolution
  /*! 
   * Nanoparticle growth one step openmp ready
  */
  int evolve_one_step_omp(double ctime);

  //! One step evolution
  /*! 
   * Nanoparticle growth one step openmp ready
  */
  int evolve_one_step_adapt(double ctime);
  
  //! Compute initial density
  /*!
   * Compute the initial density
   */
  int compute_initialdensity();
  
  //! Compute precompute
  /*! 
   * Compute constants parameters
  */
  int compute_precompute();

  //! Compute precompute
  /*! 
   * Compute constants parameters SiH4 rates
  */
  int compute_precompute_sih4();
  
  //! Compute nanoparticles potential
  /*! 
   * Compute nanoparticle potential
  */
  int compute_nanoparticle_potential();

  //! Compute collision frequencies
  /*! 
   * Compute collision frequencies
  */
  int compute_collisionfreq();

  int write_collisionfreq();

  int compute_moments(boost_array2d dens, darray &moments);
  
  int compute_moments(boost_array2d dens);
  
  int write_moments(double ctime);

  int write_partial_results(double ctime);

  int write_one_step(double ctime,
		     double begin_sim,
		     double one_step_begin);
  
  int compute_explicit_charging(double dt);

  int advance_nocharging(const double ctime);

  int advance_nocharging_omp(const double ctime);

  int advance_nocharging_ompadp(const double ctime);

  int advance_nocharging_ompadpsih4(const double ctime);

  int advance_nosplit_ompadpsih4(const double ctime);

  // fully adatptive in sections
  int advance_nosplit_sih4_adapt(const double ctime);

  // fully adatptive in sections
  int advance_nocharging_adapt(const double ctime);
  
  int compute_split_sgnucleation();

  int compute_split_sgnucleationsih4();

  int compute_nosplit_sgnucleationsih4();
  
  int compute_sgrowth_adp();
  
  int compute_sgrowth();

  int compute_sih4();
  
  int compute_sgnuc_sih4();

  int compute_coagulation();

  int compute_coagulation_omp();

  inline
  double particle_potenergy(double r, double Z) {
    return -(Kcoul*Z*eCharge*eCharge)/r;
  }

inline
double rt_affinity(double R, double Z){
    /* Computes rt r_t = \frac{Z}{\frac{Z}{R}+ A_\infty\frac{4\pi\varepsilon_0}{e^2} - \frac{5}{8}\frac{1}{R}}
    */
    double ainfinity = nm.nano.eaffinity * eCharge;
    double ai = ainfinity / (Kcoul*eCharge*eCharge);

    double rt = Z/(Z/R + ai - (5.0/(8.0*R)));
    if (rt < 0) {
        return 1000000.0;
    }
    else {
        return rt;
    }
}

  // Tunnel //FIXME
  inline
  double tunnel(double rt, double R, double Z){
      /* Griffiths Introduction to quantum mechanics The WKB approximation
          pp 295-297 eq8.24
      */
      double prefac1 = -2./GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR;
      double prefac2 = sqrt(2.*eMass*particle_potenergy(rt, Z));
      return exp(prefac1*prefac2*(rt*acos(sqrt(R/rt))-sqrt(R*(rt-R))));
  }

  inline
  double ptunnel(double Z, double R){

    double prefac1 = (-Z)*sqrt(2.0 * Kboltz
                   * cr.gm.gsys.temperature/eCharge)*(0.5/R);
    double rt = rt_affinity(R, Z);
    return prefac1*tunnel(rt, R, Z);
  }

  inline
  int omp_info() {
    BOOST_LOG_SEV(lg, info) << "Set number of threads = " << omp_get_num_threads();
    BOOST_LOG_SEV(lg, info) << "Nested = " << omp_get_nested();
    BOOST_LOG_SEV(lg, info) << "Dynamic = " << omp_get_dynamic();
    BOOST_LOG_SEV(lg, info) << "Max active levels " << omp_get_max_active_levels();
    return 0;
  }
//
// Private attributes
//
//   bgrid2d grid;                       //!< Grid vsecs x qsecs
//
//   bgrid4d grid4;                      //!< Grid vsecs x qsecs x vsecs x qsecs
//   boost_array4d efactor;              //!< Array for enhancement factor.
//   boost_array4d rcoag;                //!< Coagulation rate.
//
//   double electhermratio;          //!< Electrostatic energy to thermal ratio.
//
//   double beta0;                   //!< Coagulation rate beta0 prefactor.
//
//   //!< Container for EtaCreationFactor.
//   std::vector<EtaCreationFactor> eta_factor_vector;
//
//   //!< Container for DeathFactor.
//   std::vector<DeathFactor> death_factor_vector;
};

template < typename Type >
inline Type volume_from_radius(Type radius) {
  Type volume = 4.0*M_PI*radius*radius*radius/3.0;
  return volume;
}


template < typename Type >
inline Type radius_from_volume(Type volume) {
  return boost::math::cbrt< Type >(3.0*volume/(4.0*M_PI));
}

template < typename Type >
inline Type area_from_radius(Type radius) {
  Type area = 4.0*M_PI*radius*radius;
  return area;
}

template < typename Type >
inline Type area_from_volume(Type volume) {
  Type area = boost::math::cbrt< Type >(36.0*M_PI*volume*volume);
  return area;
}

struct qsystem{
   qsystem(){
//      NEvo n;
//      nano = n;
     }

//    qsystem(qsystem *qsys){
// //      NEvo n;
// //      nano = n;
//      }
  qsystem(const qsystem& qsys): nano(qsys.nano), l(qsys.l), tun(qsys.tun) {
    
  }
  
  qsystem(NEvo* nanoparent): nano(nanoparent) {
    tun = 0.0;
    if(nano->nm.nano.tunnel==1) {
      tun = 1.0;
    }
  }
  void operator()(const state_type &n, state_type &dndt, const value_type &t) const {

      for (unsigned int q = 1; q < nano->cr.gm.chrgs.nsections-1; ++q) {

        dndt[q] = (nano->ifreq[l][q-1]+tun*nano->tfreq[l][q-1])*n[q-1]
                + nano->efreq[l][q+1]*n[q+1]
                - n[q] *((nano->ifreq[l][q]+tun*nano->tfreq[l][q]) + nano->efreq[l][q]);
//       std::cerr << "\n dndt " << dndt[q] << "\t n " << n[q];
      }
      unsigned int q = 0;
      dndt[q] = nano->efreq[l][q+1]*n[q+1]
              - n[q] * ((nano->ifreq[l][q]+tun*nano->tfreq[l][q])); 
//       std::cerr << "\n dndt " << dndt[q] << "\t n " << n[q];

      q = nano->cr.gm.chrgs.nsections-1;
      dndt[q] = (nano->ifreq[l][q-1]+tun*nano->tfreq[l][q-1])*n[q-1]
              - n[q] * nano->efreq[l][q]; 
//       std::cerr << "\n dndt " << dndt[q] << "\t n " << n[q];
  }

  NEvo* nano;
  unsigned int l;
  double tun;
};

// inline
// int function(realtype t, N_Vector x, N_Vector dxdt, void *user_data);

// class Solver{
// public:
//   Solver(){ }

//   Solver(NEvo* nanoparent,
//          N_Vector xini_,
//          realtype ti_,
//          realtype a_=-101.,
//          realtype b_=-100.,
//          realtype reltol_=1.e-6,
//          realtype abstol_=1.e-6) :
//          nano(nanoparent),
//          xini(xini_),
//          ti(ti_),
//          a(a_),
//          b(b_),
//          reltol(reltol_),
//          abstol(abstol_) {

//     std::cerr << "\n Initial condition : "
//               << NV_Ith_S(xini,0) << '\t' << NV_Ith_S(xini,1) << '\n' << NV_LENGTH_S(xini) << "\n"; 

//     std::cerr << "\n Check : " << nano->check << "\n"; 
// //     std::cerr << "\n Initial condition : "
// //               << NV_Ith_S(xini,0) << '\t' << NV_Ith_S(xini,1) << '\n'; 
//     cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);    
//     flag = CVodeInit(cvode_mem, function, ti, xini);
//     std::cerr << "\n[ii] flag : " << flag << "\n";
//     flag = CVodeSStolerances(cvode_mem, reltol, abstol);
//     std::cerr << "\n[ii] flag : " << flag << "\n";
//     flag = CVodeSetUserData(cvode_mem, static_cast< void* >(this));
//     std::cerr << "\n[ii] flag : " << flag << "\n";

//     // ERROR WARNING check version
//     //flag = CVDense(cvode_mem, 2);
//     //std::cerr << "\n[ii] flag : " << flag << "\n";
//     //flag = CVDlsSetDenseJacFn(cvode_mem, NULL);
// //     flag = CVDiag(cvode_mem);
//     //std::cerr << "\n[ii] flag : " << flag << "\n";
//     //------------------
//     // Alternative to previous error FIXME
//     flag = CVDlsSetJacFn(cvode_mem, NULL);
//     std::cerr << "\n[ii] flag : " << flag << "\n";
//     //------------------
    
    
//     std::cerr << "\n[ii] nsecs : " << nano->cr.gm.vols.nsections << "\n";
//     // allocate vector for output
//     y = NULL;
//     y = N_VNew_Serial(2);
//   }

//   ~Solver() {
//     N_VDestroy_Serial(y);
//     CVodeFree(&cvode_mem);
//   }
//   void print() {
//     std::cerr << "\n Check : " << nano->check << "\n"; 
//   }
//   // declaring friend function grants access to a and b from class
//   friend int function(realtype t, N_Vector x, N_Vector dxdt, void *user_data);  

//   realtype a;
//   realtype b;
//   realtype reltol;
//   realtype abstol;
//   realtype ti;
//   N_Vector xini;
//   N_Vector y;

//   void compute(realtype tf, realtype dt) {
//     realtype t, tout;
//     for(iout=1, tout=ti+dt; tout <= tf; iout++, tout += dt) {
//       flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
//       std::cout << tout << '\t' << NV_Ith_S(y,0) << '\t'
//                 << NV_Ith_S(y,1) << '\t' << SUNRabs(NV_Ith_S(y,1)) << '\n'; 
//       if (flag != CV_SUCCESS) {
//         std::cerr << "\n[ee] Terminate. Flag : " << flag << "\n";
//         std::terminate();
//       }
//     }

//   }

//   void compute_step(realtype ts, realtype dt) {
//     flag = CVode(cvode_mem, ts+dt, y, &t, CV_NORMAL);
//     std::cerr << t << '\t' << NV_Ith_S(y,0) << '\t' << NV_Ith_S(y,1) << '\n'; 
//     if (flag != CV_SUCCESS) {
//       std::cerr << "\n[ee] Terminate. Flag : " << flag << "\n";
//       std::terminate();
//     }
//   }

// private:
//   void *cvode_mem = NULL;
//   int flag;
//   realtype t;
//   int iout;
//   NEvo *nano;
// };

// // function to integrate
// int function(realtype t, N_Vector x, N_Vector dxdt, void *user_data) {

//   // static_cast user_data to class Solver 
//   Solver* s = static_cast< Solver* >(user_data); 

//   realtype x0, x1;

//   x0 = NV_Ith_S(x, 0);
//   x1 = NV_Ith_S(x, 1);

//   NV_Ith_S(dxdt, 0) = s->a*x0 + s->b*x1; // in boost dxdt[ 0 ] = -101.0 * x[ 0 ] - 100.0 * x[ 1 ];
//   NV_Ith_S(dxdt, 1) = x0; // in boost dxdt[ 1 ] = x[ 0 ];

//   return(0);
// }

// ------------------------------ lsoda solver

inline
int ls_qsystemf(double t, double *n, double *dndt, void *data);

struct ls_qsystem{
  ls_qsystem() {}

  ls_qsystem(const qsystem& qsys): nano(qsys.nano), l(qsys.l), tun(qsys.tun) {}
  
  ls_qsystem(NEvo* nanoparent): nano(nanoparent) {
    tun = 0.0;
    if(nano->nm.nano.tunnel==1) {
      tun = 1.0;
    }
  }

  friend int ls_qsystemf(double t, double *n, double *dndt, void *data);

  NEvo* nano;
  unsigned int l;
  double tun;
};

int ls_qsystemf(double t, double *n, double *dndt, void *data) {

  const ls_qsystem* qsys = static_cast< ls_qsystem* >(data);

  unsigned int nchrgs = qsys->nano->cr.gm.chrgs.nsections;
  unsigned int l = qsys->l;
  
  for (unsigned int q = 1; q < nchrgs-1; ++q) {

    dndt[q] = (qsys->nano->ifreq[l][q-1]+qsys->tun*qsys->nano->tfreq[l][q-1])*n[q-1]
      + qsys->nano->efreq[l][q+1]*n[q+1]
      - n[q] *((qsys->nano->ifreq[l][q]+qsys->tun*qsys->nano->tfreq[l][q]) + qsys->nano->efreq[l][q]);

  }
  unsigned int q = 0;
  dndt[q] = qsys->nano->efreq[l][q+1]*n[q+1]
    - n[q] * ((qsys->nano->ifreq[l][q]+qsys->tun*qsys->nano->tfreq[l][q])); 

  q = nchrgs-1;
  dndt[q] = (qsys->nano->ifreq[l][q-1]+qsys->tun*qsys->nano->tfreq[l][q-1])*n[q-1]
    - n[q] * qsys->nano->efreq[l][q];

  return 0;
}

// ------------------------------ lsoda solver


inline
double gaussian_distribution(double x, double mu, double sigma){
  // mu -> mean
  // sigma -> standard deviation
  double prefac = 1.0/sqrt(2.0*pi*sigma*sigma);
  double arg = -0.5*pow((x-mu)/sigma, 2);
  return prefac*exp(arg);
}

#endif // NEVO_H
