/*
 * Copyright 2017 Benjamin Santos <caos21@gmail.com>
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

#ifndef EINT_H
#define EINT_H

#define FINITE 1

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>

#include <omp.h>

#include <boost/array.hpp>
#include <boost/math/tools/roots.hpp>
// #include <boost/numeric/odeint.hpp>

//#include <boost/math/special_functions/factorials.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_sparse.hpp>
#include <boost/numeric/bindings/lapack/gesv.hpp>
#include <boost/numeric/bindings/lapack/sysv.hpp>
#include <boost/numeric/bindings/traits/ublas_vector2.hpp>

// #include <boost/numeric/bindings/umfpack/umfpack.hpp>

#include "constants.h"

namespace tools = boost::math::tools;
namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;
// namespace umf = boost::numeric::bindings::umfpack;
namespace bmath = boost::math;

#include <boost/multiprecision/cpp_int.hpp>
using boost::multiprecision::cpp_int;
#include <boost/multiprecision/cpp_dec_float.hpp>
using boost::multiprecision::cpp_dec_float_50;
using boost::multiprecision::cpp_dec_float_100;

typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<50>> mpfloat;

//typedef double mpfloat;

//typedef long double mpfloat;

struct TerminationCondition  {
  bool operator() (double min, double max)  {
    return abs((min - max)/min) <= 0.0001;
  }
};  

namespace eint {
  typedef ublas::matrix<double> ubmatrix;
  typedef ublas::vector<double> ubvector;

  //! Factorized (dimensionless) potential approximation IPA.
  /*!
    \param rt the separation.
    \param r21 the ratio r2/r1.
    \param q21 the ratio q2/q1.
    \param eps the dielectric constant.
  */
  inline
  double potential_ipa_fact(const double rt, const double r21,
                            const double q21, const double eps) {
    double kappa = (eps-1.0)/(eps+2.0);
    double A = kappa*Kcoul*pow(r21, 3);
    double B = kappa*Kcoul*q21*q21;
    return Kcoul*q21/rt
          - A/(2*rt*rt*(rt*rt - r21*r21))
          - B/(2*rt*rt*(rt*rt - 1.0));
  }

  //! Factorized (dimensionless) force approximation IPA.
  /*!
    \param rt the separation.
    \param r21 the ratio r2/r1.
    \param q21 the ratio q2/q1.
    \param eps the dielectric constant.
  */
  inline
  double force_ipa_fact(const double rt, const double r21,
                        const double q21, const double eps) {
    double kappa = (eps-1.0)/(eps+2.0);
    double A = kappa*Kcoul*pow(r21, 3) * (2.*rt*rt-r21*r21);
    double B = kappa*Kcoul*q21*q21*(2.*rt*rt-1.);
    return Kcoul*q21/(rt*rt)
          - A/(pow(rt, 3)*pow((rt*rt - r21*r21),2))
          - B/(pow(rt, 3)*pow(rt*rt - 1.0, 2));
  }


  /************/ 
  //! Series potential approximation.
  /*!
  * Potential energy for a point charge q and dielectric sphere radius a
  * at a distance s, Stratton pag. 204
  * Here we add the interaction of the other particle. Superposed.
  * SIPA
    \param rt the separation.
    \param r21 the ratio r2/r1.
    \param q21 the ratio q2/q1.
    \param eps the dielectric constant.
    \param terms the number of terms in the expansion.
  */
  inline
  double potential_series(const double r1, const double q1,
                          const double r2, const double q2, const double rt,
                          const double eps, unsigned int terms=50) {

      // Potential r1<<r2 q2=0
      double prefactor1 = Kcoul * q1*q1 * (eps-1.0) / (2.*r2);
      double pot1 = 0.0;
      for(unsigned int j=0.0; j<terms; ++j) {
        double numer = j*pow(r2/rt, 2*j+2);
        double denom = (j*eps+j+1.0);
        pot1 += numer/denom;
      }
      pot1 *= prefactor1;

      //Potential r2<<r1 q1=0
      double prefactor2 = Kcoul * q2*q2 * (eps-1.0) / (2.*r1);
      double pot2 = 0.0;
      for(unsigned int j=0.0; j<terms; ++j) {
        double numer = j*pow(r1/rt, 2*j+2);
        double denom = (j*eps+j+1.0);
        pot2 += numer/denom;
      }
      pot2 *= prefactor2;

      //Coulomb potential energy
      double pot_c = Kcoul*q1*q2/rt;

      return pot_c - pot1 - pot2;
  }

  //! Series force approximation.
  /*!
  * Force for a point charge q and dielectric sphere radius a
  * at a distance s, Stratton pag. 204
  * Here we add the interaction of the other particle. Superposed.
  * SIPA
    \param rt the separation.
    \param r21 the ratio r2/r1.
    \param q21 the ratio q2/q1.
    \param eps the dielectric constant.
    \param terms the number of terms in the expansion.
  */
  inline
  double force_series(const double r1, const double q1,
                          const double r2, const double q2, const double rt,
                          const double eps, unsigned int terms=50) {

      // Potential r1<<r2 q2=0
      double prefactor1 = Kcoul * q1*q1;
      double force1 = 0.0;
      for(unsigned int j=0.0; j<terms; ++j) {
        double numer = pow(r2, 2*j+1)*(j*(j+1)*(eps-1.0));
        double denom = pow(rt, 2*j+3)*(j*eps+j+1.0);
        force1 += numer/denom;
      }
      force1 *= prefactor1;

      //Potential r2<<r1 q1=0
      double prefactor2 = Kcoul * q2*q2;
      double force2 = 0.0;
      for(unsigned int j=0.0; j<terms; ++j) {
        double numer = pow(r1, 2*j+1)*(j*(j+1)*(eps-1.0));
        double denom = pow(rt, 2*j+3)*(j*eps+j+1.0);
        force2 += numer/denom;
      }
      force2 *= prefactor2;

      //Coulomb potential energy
      double force_c = Kcoul*q1*q2/(rt*rt);

      return force_c - force1 - force2;
  }
  /************/

  // unsigned int factorial(unsigned int n) {
  //   return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
  // }
  inline
  double factorial(unsigned int uin) {
    double n = static_cast<double>(uin);
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
  }

  template <typename T>
  T from_stirling(double x){
    // helper, derived from Stirling's approximation
    T y = x;
    T ret = y*log(y)+0.5*log(2.*pi*y);
    return ret;
  }

  template <typename T>
  T stirling_quotient(double x, double y){
    T ret = from_stirling<T>(x+y) - from_stirling<T>(x) - from_stirling<T>(y);
    return ret;
  }

  // compute D matrix elements
  inline
  double dijlm(double rt, double r1, double r2, unsigned int l, unsigned int m) {
    if (m==0) {
      return pow(r1/rt, l);
    }
    
    unsigned int lm = l+m;
    if (lm <= NFMAX-1){
      // use exact factorials
      return (Facts[lm]/(Facts[l]*Facts[m]))*pow(r1/rt, l)*pow(r2/rt, m);
    }
    else {
      // Stirling's approximation
      double lpart = l*(log((lm)*r1/(l*rt)));
      double mpart = m*(log((lm)*r2/(m*rt)));
      double apart = 0.5*log(lm/(2.*pi*l*m));
      double hoterms = 1./(12*lm) - 1./(360*lm*lm*lm);
      return exp(lpart+mpart+apart+hoterms);
    }
  }

  // compute C Matrix
  inline
  double clm(double rt, double r21, unsigned int l, unsigned int m, double eps) {
    if (l==m) {
      double cval = (rt/r21)*(((eps+1)*l+1)/((eps-1)*l));
      return cval;
    }
    return 0.0;
  }
  
  // compute multipolar coefficients
  inline
  void compute_MPCoefficients(ubvector &A1coefficients,
			      ubvector &A2coefficients,
			      double &A10,
			      double &A20,
			      const double rt,
			      const double r21,
			      const double q21,
			      const double eps,
			      const unsigned int nterms) {

    A10 = Kcoul;
    A20 = Kcoul*q21/r21;

    // system is nterms - 0th term
    unsigned int ssize = nterms-1;

#ifdef USING_EIGEN

    VectorDD b1 = VectorDD::Zero(ssize);
    VectorDD b2 = VectorDD::Zero(ssize);

    MatrixDD C1 = MatrixDD::Zero(ssize, ssize);
    MatrixDD C2 = MatrixDD::Zero(ssize, ssize);

    MatrixDD C1_inv = MatrixDD::Zero(ssize, ssize);
    MatrixDD C2_inv = MatrixDD::Zero(ssize, ssize);
    
    MatrixDD D1 = MatrixDD::Zero(ssize, ssize);
    MatrixDD D2 = MatrixDD::Zero(ssize, ssize); 
    
#else
    ubvector b1 = ublas::zero_vector<double>(ssize);
    ubvector b2 = ublas::zero_vector<double>(ssize);

    ubmatrix C1 = ublas::zero_matrix<double>(ssize, ssize);
    ubmatrix C2 = ublas::zero_matrix<double>(ssize, ssize);

    ubmatrix C1_inv = ublas::zero_matrix<double>(ssize, ssize);
    ubmatrix C2_inv = ublas::zero_matrix<double>(ssize, ssize);
    
    ubmatrix D1 = ublas::zero_matrix<double>(ssize, ssize);
    ubmatrix D2 = ublas::zero_matrix<double>(ssize, ssize);
#endif
    
    for (unsigned int l=0; l<ssize; ++l) {
      // b vectors
      b1(l) = -dijlm(rt, 1.0, r21, l+1, 0) * A20;
      b2(l) = -dijlm(rt, r21, 1.0, l+1, 0) * A10;
      
      //Fill C Matrices (diagonal)
      C1(l, l) = clm(rt, r21, l+1, l+1, eps);
      C2(l, l) = clm(rt, 1.0, l+1, l+1, eps);

      //inverse C diagonal matrices
      C1_inv(l, l) =  1.0/C1(l, l);
      C2_inv(l, l) =  1.0/C2(l, l);
    }
    for (unsigned int l=0; l<ssize; ++l) {
      for (unsigned int m=0; m<ssize; ++m) {
	//Fill D Matrices
	D1(l, m) = dijlm(rt, 1.0, r21, l+1, m+1);
	D2(l, m) = dijlm(rt, r21, 1.0, l+1, m+1);
      }
    }

#ifdef USING_EIGEN
 #ifdef USING_AMPC//using complete matrix A
    MatrixDD M = MatrixDD::Zero(2*ssize, 2*ssize);

    M.block(0    ,       0, ssize  ,   ssize) = C1;
    M.block(0    ,   ssize, ssize  ,   ssize) = D1;
    
    M.block(ssize,       0, ssize  ,   ssize) = C2;
    M.block(ssize,   ssize, ssize  ,   ssize) = D2;

    VectorDD B = VectorDD::Zero(2*ssize);

    B.segment(0    , ssize) = b1;
    B.segment(ssize, ssize) = b2;

    /* MatrixDD S = M.selfadjointView<Eigen::Upper>(); */

    /* Eigen::ConjugateGradient<MatrixDD, Eigen::Upper> cg; */
    /* cg.compute(S); */

    /* VectorDD X(2*ssize); */

    /* X = cg.solve(B); */
    
    /* std::cerr << "[ii] #iterations:     " << cg.iterations() << std::endl; */
    /* std::cerr << "[ii] estimated error: " << cg.error()      << std::endl; */

    VectorDD X = VectorDD::Zero(2*ssize);
    X = M.colPivHouseholderQr().solve(B);

    for (unsigned int l=0; l<ssize; ++l) {
      A1coefficients(l) = X(l);
      A2coefficients(l) = X(l+ssize);
    }
 #else//USING_AMPC
    MatrixDD M12 = MatrixDD::Zero(ssize, ssize);
    VectorDD B1 = VectorDD::Zero(ssize);

    M12 = C1 - D1 * (C2_inv * D2);
    B1 =  b1 - D1 * (C2_inv * b2);

    VectorDD X1 = VectorDD::Zero(ssize);
    //X1 = M12.colPivHouseholderQr().solve(B1);
    X1 = M12.partialPivLu().solve(B1);

    /* MatrixDD M21 = MatrixDD::Zero(ssize, ssize); */
    /* VectorDD B2 = VectorDD::Zero(ssize); */

    /* M21 = C2 - D2 * (C1_inv * D1); */
    /* B2 =  b2 - D2 * (C1_inv * b1); */

    VectorDD X2 = VectorDD::Zero(ssize);
    //X2 = M21.colPivHouseholderQr().solve(B2);
    // X2 = C2inv * (b2 - D2 * X1)
    X2 = C2_inv * (b2 - D2 * X1);
    
    for (unsigned int l=0; l<ssize; ++l) {
      A1coefficients(l) = X1(l);
      A2coefficients(l) = X2(l);
    }
 #endif//USING_AMPC
#else//USING_EIGEN
 #ifdef USING_AMPC//using complete matrix A
    //
    // ************  Slower, using full matrix ****************
    //
    ubmatrix M = ublas::zero_matrix<double>(2*ssize, 2*ssize);

    ublas::range r0s = ublas::range(0, ssize);
    ublas::range rs2s = ublas::range(ssize, 2*ssize);

    project(M, r0s, r0s) = C1;
    project(M, r0s, rs2s) = D1;

    project(M, rs2s, r0s) = C2;
    project(M, rs2s, rs2s) = D2;

    ubvector B = ublas::zero_vector<double>(2*ssize);

    project(B, r0s) = b1;
    project(B, rs2s) = b2;

    ublas::symmetric_adaptor<ublas::matrix<double>, ublas::lower> SM(M);
    lapack::sysv(SM, B);
    //lapack::gesv(M, B);

    A1coefficients = project(B, r0s);
    A2coefficients = project(B, rs2s);

 #else//USING_AMPC using complete matrix A
    //
    // ************  Using reduced system  ****************
    //    
    ubmatrix M12 = ublas::zero_matrix<double>(ssize, ssize);
    ubvector B1 = ublas::zero_vector<double>(ssize);

    // WARNING future optimisation prod(A, temp_type(prod(B,C));

    //M12 = C1 - D1 * ublas::prod(C2_inv, D2);
    ubmatrix MTemp = ublas::prod(C2_inv, D2);
    M12 = C1 - ublas::prod(D1, MTemp);

    // B1 =  b1 - D1 * (C2_inv * b2);
    ubvector VTemp = ublas::prod(C2_inv, b2);
    B1 =  b1 - ublas::prod(D1, VTemp);
    
    ubmatrix M21 = ublas::zero_matrix<double>(ssize, ssize);
    ubvector B2 = ublas::zero_vector<double>(ssize);

    //M21 = C2 - D2 * ublas::prod(C1_inv, D1);
    MTemp = ublas::prod(C1_inv, D1);
    M21 = C2 - ublas::prod(D2, MTemp);

    //B2 =  b2 - D2 * (C1_inv * b1);
    VTemp = ublas::prod(C1_inv, b1);
    B2 =  b2 - ublas::prod(D2, VTemp);

    //    lapack::gesv(M12, B1);
    // solve symmetric linear system of equations
    /* lapack::gesv(M12, B1); */
    ublas::symmetric_adaptor<ublas::matrix<double>, ublas::lower> SM12(M12);
    lapack::sysv(SM12, B1);

    // results
    A1coefficients = B1;

    // solve symmetric linear system of equations
    ublas::symmetric_adaptor<ublas::matrix<double>, ublas::lower> SM21(M21);
    lapack::sysv(SM21, B2);

    A2coefficients = B2;    
    /* ubvector VTemp2  = ublas::zero_vector<double>(ssize); // = ublas::prod(D2, A1coefficients); */
    /* ublas::axpy_prod(D2, A1coefficients, VTemp2, true); */
    /* VTemp = b2 - VTemp2; */
    /* //A2coefficients = ublas::prod(C2_inv, VTemp); */
    /* ublas::axpy_prod(C2_inv, VTemp, A2coefficients, true); */
    
    
    //    std::cerr << "\n A1 " << A1coefficients << '\n';
    /* std::cerr << "\n A2 " << A2coefficients << '\n'; */
 #endif//USING_AMPC
#endif//USING_EIGEN
  }

  //! Bichoutskaia potential.
  /*!
   *
   *As defined in equation A7 of
    Lindgren, E. B., Chan, H.-K., Stace, A. J. & Besley, E.
    Progress in the theory of electrostatic interactions between charged particles.
    Phys. Chem. Chem. Phys. 18, 5883–5895 (2016)
    \param r1 radius of particle 1.
    \param r2 radius of particle 2.
    \param q1 charge of particle 1.
    \param q2 charge of particle 2.
    \param h separation between particles.
    \param eps the dielectric constant.
    \param acoefficients the multipole moment coefficients
    \
  */
  inline
  double mpc_potential(const double& rt, const double& r21, const double& q21,
                       const ubvector& A1coefficients, const ubvector& A2coefficients,
		       const double& eps) {

    unsigned int size = A1coefficients.size();

    double invK = 1.0/Kcoul;
    
    double potc = Kcoul*q21/rt;
    
    double epsm = eps-1.0;
    double epsp = eps+1.0;
    
    double pot1 = 0.0;
    for(unsigned int mp=1; mp<size+1; ++mp) {
      pot1 += A2coefficients(mp-1)*pow(r21/rt, mp+1);
      //std::cerr << "\n pot1 " << pot1;
    }
    //std::cerr << "\n";

    double pot2 = 0.0;
    for(unsigned int lp=1; lp<size+1; ++lp) {
      double prefac = (epsp*lp+1)/(epsm*lp);
      pot2 += prefac * A1coefficients(lp-1) * A1coefficients(lp-1);
    }    
    // NOTE 1/2 factor correction to potential as stated in
    // 1.Stace, A. J. & Bichoutskaia, E. Reply to the ‘Comment on “Treating
    // highly charged carbon and fullerene clusters as dielectric particles”’
    // by H. Zettergren and H. Cederquist, Phys. Chem. Chem. Phys., 2012, 14,
    // DOI: 10.1039/c2cp42883k. Phys. Chem. Chem. Phys. 14, 16771–16772 (2012).
    /* std::cerr << '\n' << pot_coul << '\t' << pot_2 << '\t' << pot_3 << '\n'; */
    return potc + 0.5*pot1 - 0.5*invK*pot2;
  }

  
  //! Bichoutskaia potential.
  /*!
   *
   *As defined in equation A7 of
    Lindgren, E. B., Chan, H.-K., Stace, A. J. & Besley, E.
    Progress in the theory of electrostatic interactions between charged particles.
    Phys. Chem. Chem. Phys. 18, 5883–5895 (2016)
    \param r1 radius of particle 1.
    \param r2 radius of particle 2.
    \param q1 charge of particle 1.
    \param q2 charge of particle 2.
    \param h separation between particles.
    \param eps the dielectric constant.
    \param acoefficients the multipole moment coefficients
    \
  */
  inline
  double potential_bichoutskaia(const double r1, const double q1,
                                const double r2, const double q2,
                                const double h,
                                const ubvector& acoefficients,
                                const double eps) {

    // WARNING
    return 1.0;
  }

  //! Bichoutskaia force.
  /*!
  * 
  *As defined in equation 5 of 
    Lindgren, E. B., Chan, H.-K., Stace, A. J. & Besley, E.
    Progress in the theory of electrostatic interactions between charged particles.
    Phys. Chem. Chem. Phys. 18, 5883–5895 (2016)

    This is the scaled force.

    \param A10 Multipolar coefficent 0.
    \param A1coefficients Multipolar coefficients
  */
  inline
  double mpc_force(const double& A10,
		   const ubvector& A1coefficients,
		   double eps) {
    
    unsigned int size = A1coefficients.size();
    
    double invK = 1.0/Kcoul;
    
    double epsm = eps-1.0;
    double epsp = eps+1.0;
    
    double force = 0.0;

    // Populate new vector with A1 coefficients
    ubvector Acoeffs = ublas::zero_vector<double>(size+1);
    
    for(unsigned int m=1; m<size+1; ++m) {
      Acoeffs(m) = A1coefficients(m-1);
    }
    Acoeffs(0) = A10;
    
    for(unsigned int l=0; l<size; ++l) {
      double prefac = (epsp*(l+1)+1)/epsm;
      force += prefac*Acoeffs(l)*Acoeffs(l+1);
    }
    
    return -invK * force;
  }

  //! Bichoutskaia force.
  /*!
  * \param r1 radius of particle 1.
    \param r2 radius of particle 2.
    \param q1 charge of particle 1.
    \param q2 charge of particle 2.
    \param h separation between particles.
    \param eps the dielectric constant.
    \param acoefficients the multipole moment coefficients
    \
  */
  inline
  double force_bichoutskaia(const double r1, const double q1,
                            const double r2, const double q2,
                            const double h,
                            const ubvector& acoefficients,
                            const double eps) {
    // WARNING to complete
    return -1.0;
  }
  
  inline
  double relative_error(const double trueval, const double expval){
    return abs((expval-trueval)/trueval);
  }

  inline
  double max_relative_error(const double a, const double b){
    return abs((a-b)/std::min(a,b));
  }
  
  struct potential_ipa_funct
  {
    potential_ipa_funct(double r21_, double q21_, double eps_):
                        r21(r21_), q21(q21_), eps(eps_){
    }

    double operator()(double const& rt) {
      return potential_ipa_fact(rt, r21, q21, eps);
    }

    double r21;
    double q21;
    double eps;
  };

  struct force_ipa_funct
  {
    force_ipa_funct(double r21_, double q21_, double eps_):
                        r21(r21_), q21(q21_), eps(eps_){
    }

    double operator()(double const& rt) {
      return force_ipa_fact(rt, r21, q21, eps);
    }

    double r21;
    double q21;
    double eps;
  };




  // multipolar coefficients potential functor
  struct potential_mpc_funct
  {
    potential_mpc_funct(double r21_,
                        double q21_,
                        double eps_,
                        unsigned int nterms_):
                        r21(r21_), q21(q21_), eps(eps_),
                        nterms(nterms_) {

      // system is nterms - 0th term
      unsigned int ssize = nterms-1;
    
      // 1, n-1 terms
      A1coefficients = ublas::zero_vector<double>(ssize);
      A2coefficients = ublas::zero_vector<double>(ssize);

      // 0 terms
      A10 = 0.0;
      A20 = 0.0;
    }

    // constructor if coefficients are known
    potential_mpc_funct(const ubvector &A1coefficients_,
			const ubvector &A2coefficients_,
			double A10_,
			double A20_,
			double eps_,
			unsigned int nterms_):
                        A1coefficients(A1coefficients_),
		        A2coefficients(A2coefficients_),
		        eps(eps_),
                        nterms(nterms_) {

    }
    
    double operator()(double const& rt) {

      // compute coefficients
      compute_MPCoefficients(A1coefficients, A2coefficients,
			     A10, A20, rt, r21, q21, eps, nterms);

      // compute Bichoutskaia force
      return mpc_potential(rt, r21, q21, A1coefficients, A2coefficients, eps);
    }

    ubvector A1coefficients;
    ubvector A2coefficients;
    double A10;
    double A20;
    
    double r21;
    double q21;
    double eps;
    unsigned int nterms;   

  };

  // multipolar coefficients force functor
  struct force_mpc_funct {
    force_mpc_funct(double r21_,
                    double q21_,
                    double eps_,
                    unsigned int nterms_):
                    r21(r21_), q21(q21_), eps(eps_),
                    nterms(nterms_) {

      // system is nterms - 0th term
      unsigned int ssize = nterms-1;

      // 1, n-1 terms
      A1coefficients = ublas::zero_vector<double>(ssize);
      A2coefficients = ublas::zero_vector<double>(ssize);

      // 0 terms
      A10 = 0.0;
      A20 = 0.0;
    }

    // constructor if coefficients are known
    force_mpc_funct(const ubvector &A1coefficients_,
		    const ubvector &A2coefficients_,
		    double A10_,
		    double A20_,
                    double eps_,
                    unsigned int nterms_):
                      A1coefficients(A1coefficients_),
		      A2coefficients(A2coefficients_),			
		      eps(eps_),
                      nterms(nterms_) {

    }
    
    double operator()(double const& rt) {

      // compute coefficients
      compute_MPCoefficients(A1coefficients, A2coefficients,
      			     A10, A20, rt, r21, q21, eps, nterms);

      // compute Bichoutskaia force
      return mpc_force(A10, A1coefficients, eps);

    }


    ubvector A1coefficients;
    ubvector A2coefficients;
    double A10;
    double A20;
    
    double r21;
    double q21;
    double eps;
    unsigned int nterms;

#ifdef USING_EIGEN

    /* MatrixDD  */
    
#endif

  };
  
  // compute eta factor for a pair of particles
  inline
  double efactor_mpc(double r1, double r2,
                     double q1, double q2,
                     const double eps, const double temperature,
                     const unsigned int nterms) {

    // for numerical stability, accuracy and fast convergence we want r21 < 1.0
    if(r2>r1){
      double raux = r2;
      r2 = r1;
      r1 = raux;
      double qaux = q2;
      q2 = q1;
      q1 = qaux;
    }

    double r21 = r2/r1;

    // WARNING float comparison
    // avoid division by 0
    if(q1==0.0){
      // q1=0, permute particle 1 with 2
      double raux = r2;
      r2 = r1;
      r1 = raux;
      double qaux = q2;
      q2 = q1;
      q1 = qaux;
      r21 = r2/r1;
    }

    double q21 = q2/q1;

    // max rx
    double max = 1.0e-5/r1;

    // scaled contact radii
    double rt = 1.0 + r21;
    double min = rt;

    // Find scaled potential at contact and nterms
    unsigned int curr_nterms = nterms;
    unsigned int step_nterms = 10;
    unsigned int iter = 0;
    unsigned int max_iter = 20;

    double pmpc_rt;
    double pmpc_rt_prev = 0.0;
    
    while (true) {      
      potential_mpc_funct pmpcfunct(r21, q21, eps,
				    curr_nterms);

      pmpc_rt = pmpcfunct(rt);

      if ((iter>0) && (max_relative_error(pmpc_rt, pmpc_rt_prev)<PotRE)) {
	std::cerr << "\n[II] Convergence achieved for N = " << curr_nterms;
	break;
      }
      if (iter>max_iter) {
	std::cerr << "\n[EE] Max iterations = " <<  iter;
	std::cerr << "\n[EE] Number of terms = " << curr_nterms;
	break;
      }
      
      pmpc_rt_prev = pmpc_rt;
      curr_nterms += step_nterms;
      ++iter;
    }
    

    // WARNING copy instance
    potential_mpc_funct pmpcfunct(r21, q21, eps,
				  curr_nterms);
    
    double potprefactor = q1*q1*eCharge*eCharge/r1;

    // Potential at contact
    double phi_rt = potprefactor*pmpc_rt;

    if(pmpc_rt>0){
      double eta = exp(-phi_rt/(Kboltz*temperature));
      std::cerr << "\n[ii] Repulsive phi = " << phi_rt;
      /* if(phi_rt < 0){ */
      /* 	std::terminate(); */
      /* } */
      return eta;
    }
    else {// potential is negative (attractive) or zero at contact

      // force functor
      force_mpc_funct forcempcfunct(r21, q21, eps, curr_nterms);
      
      // Force at contact
      double forcempc_rt = forcempcfunct(rt);
      
      // Force at r max
      double forcempc_max = forcempcfunct(max);
      
      // checks if force is monotonically decreasing [non bracketed]
      if(forcempc_rt*forcempc_max >= 0.0){
	std::cerr << "\n[ii] Attractive phi = " << phi_rt;
	double eta = 1.0 - phi_rt /(Kboltz*temperature);
	/* if(phi_rt > 0){ */
	/*   std::terminate(); */
	/* } */
	return eta;
      }
      else {// find maximum, zero of force is bracketed between rt and max
	std::cerr << "\n[ii] Mixed phi_rt = " << phi_rt;
	std::pair<double, double> pair_pmpc;
	boost::uintmax_t bmax_iter = 1000;
	tools::eps_tolerance<double> tol(30);
	try {
	  //#pragma omp critical
	  //{	  
	  pair_pmpc = tools::toms748_solve(forcempcfunct, min, max, forcempc_rt, forcempc_max, tol, bmax_iter);

	  //pair_pmpc = tools::toms748_solve(forcempcfunct, rt, max, tol, max_iter);
	
	 // std::cerr << "\n" << rt << "\t" << max << "\t" << forcempcfunct(rt) << "\t" << forcempcfunct(max);
	 //}
	}
	catch(const std::exception& exc) {
	  
	  //pair_pmpc.first = 0.0;
	  //pair_pmpc.second = 0.0;

	    std::cerr << '\n' << exc.what() << '\n';

	  //int nst = 50;
	  //double st = (max-rt)/(nst-1);
	  //for(int i=0; i<50; ++i){
	  //  double rr = rt + st*i;
	  //  std::cerr << "\n" << rr << "\t" << pmpcfunct(rr) << "\t" << forcempcfunct(rr);
	  //}
	  //std::cerr << '\n';	
	    std::terminate();
	  //return 1.0 - phi_rt /(Kboltz*temperature);
	}
	if(max_iter > 990){
	  std::cerr << "\n ERROR max iter " << max_iter << "\n\n";
	  std::terminate();
	}
	double rbarrier = 0.5*(pair_pmpc.first+pair_pmpc.second);
	if(rbarrier>=0){
	  double phimax = potprefactor * pmpcfunct(rbarrier);
	  double eta = exp(-phimax/(Kboltz*temperature))
	               *(1.0+(phimax-phi_rt)/(Kboltz*temperature));	  
	  return eta;
	}
	else{
	  std::cerr << "\n ERROR Negative rbarrier " << rbarrier << '\n';
	  std::terminate();
	}
      }
    }
    
    return 1.0;
  }


  // compute ipa enhancement factor for a pair of particles
  // WARNING using pmpc instead of pipa
  inline
  double efactor_ipa2(double r1, double r2,
                      double q1, double q2,
                      double eps, double temperature) {

    double r21 = r2/r1;
    double q21 = q2/q1;
    // check q1==0
    if(fabs(q1)<1.e-200){
      // q1=0, permute particle 1 with 2
      double raux = r2;
      r2 = r1;
      r1 = raux;
      double qaux = q2;
      q2 = q1;
      q1 = qaux;
      r21 = r2/r1;
      q21 = 0.0;
    }

    double max = 500.0;
    boost::uintmax_t max_iter = 1000;
    tools::eps_tolerance<double> tol(30);

    double rt = 1.0 + r21;

    double min = rt;

    // potential function
    potential_ipa_funct pipafunct(r21, q21, eps);
    //double pipa = pipafunct(rt);

    potential_mpc_funct pmpcfunct(r21, q21, eps, 25);
    double pmpc = pmpcfunct(rt);

    // force function
    force_ipa_funct forceipafunct(r21, q21, eps);

    std::pair<double, double> pair_pipa;
    bool failed = false;
    try {
//       pair_pipa = tools::toms748_solve(pipafunct, min, max, tol, max_iter);
      pair_pipa = tools::toms748_solve(forceipafunct, min, max, tol, max_iter);
    }
    catch(const std::exception& e) {
      failed = true;
      pair_pipa.first = 0.0;
//       std::cout << "\nMessage from thrown exception was:\n  " << e.what() << std::endl;
//       std::cout << "\nmax_iter = " << max_iter;
    }
    bool attcontact = (pmpc<0.0? true: false);
    bool fullatt = attcontact && failed;
    bool withphimax = !failed || attcontact;

    double potprefactor = q1*q1*eCharge*eCharge/r1;
    double eta = 0.0;
    double phimin = potprefactor*pmpc;

    if(withphimax){
//       double phimax = potprefactor*pair_pipa.first;
      double phimax = 0.0;
      if (pair_pipa.first>0.0){
        // NOTE could be better midpoint
        // pair_pipa.first + (pair_pipa.second - pair_pipa.first)/2;
        phimax = potprefactor*pmpcfunct(pair_pipa.first);
      }

      eta = exp(-phimax/(Kboltz*temperature))
          *(1.0+(phimax-phimin)/(Kboltz*temperature));
//       std::cout << std::endl
//               << "\t pipa.first = " << pair_pipa.first
//               << "\t pipa.second = " << pair_pipa.second
//               << "\t max_iter = " << max_iter;

    }
    if(fullatt){
      eta = 1.0 - phimin /(Kboltz*temperature);
    }
    if(!attcontact){
      eta = exp(-phimin/(Kboltz*temperature));
    }

//     std::cout << std::endl << "r21 = " << r21 << "\tq21 = " << q21 
//               << "\t phimax = " << potprefactor*pair_pipa.first
//               << "\t phimin = " << phimin << "\tcontact = " << attcontact
//               << "\tfullatt = " << fullatt << "\teta = " << eta;
    return eta;
  }

  // compute ipa enhancement factor for a pair of particles
  // WARNING using pmpc instead of pipa
  // we use ipa force to find the location of the barrier if
  // the mpc potential is attractive at contact
  // then we compute the enhancement factor using the mpc potential
  // at contact and at the barrier determined by zeroing the ipa force
  inline
  double efactor_ipa(double r1, double r2,
                     double q1, double q2,
                     const double eps, const double temperature,
                     const unsigned int nterms) {


    if(r2>10.0*r1){
      double raux = r2;
      r2 = r1;
      r1 = raux;
      double qaux = q2;
      q2 = q1;
      q1 = qaux;
    }

    
    double r21 = r2/r1;
    double q21 = q2/q1;

    // WARNING float comparison
    if(q1==0.0){
      // q1=0, permute particle 1 with 2
      double raux = r2;
      r2 = r1;
      r1 = raux;
      double qaux = q2;
      q2 = q1;
      q1 = qaux;
      r21 = r2/r1;
      q21 = 0.0;
    }
    
    
    double max = 1.0e-5/r1;

    double rt = 1.0 + r21;

    double min = rt;

    // // ipa functor
    // potential_ipa_funct pipafunct(r21, q21, eps);

    // // ipa at contact 
    // double pipa_rt = pipafunct(rt)

    // mpc functor
    potential_mpc_funct pmpcfunct(r21, q21, eps, nterms);

    // mpc at contact
    double pmpc_rt = pmpcfunct(rt);

    // potential prefactor
    double potprefactor = q1*q1*eCharge*eCharge/r1;

    // Potential at contact in mks
    double phi_rt = potprefactor*pmpc_rt;

    // repulsive potential
    if(pmpc_rt>0){
      double eta = exp(-phi_rt/(Kboltz*temperature));
      /* if(phi_rt < 0){ */
      /* 	std::terminate(); */
      /* } */
      return eta;
    }
    else {// potential is negative (attractive) or zero at contact

      // force functor
      force_ipa_funct forceipafunct(r21, q21, eps);
      
      // Force at contact
      double forceipa_rt = forceipafunct(rt);
      
      // Force at r max
      double forceipa_max = forceipafunct(max);
      
      // checks if force is monotonically decreasing [non bracketed]
      if(forceipa_rt*forceipa_max >= 0.0){
	double eta = 1.0 - phi_rt /(Kboltz*temperature);
	/* if(phi_rt > 0){ */
	/*   std::terminate(); */
	/* } */
	return eta;
      }
      else {// find maximum, zero of force is bracketed between rt and max
	std::pair<double, double> pair_pipa;
	boost::uintmax_t max_iter = 1000;
	tools::eps_tolerance<double> tol(30);
	try {
	  //#pragma omp critical
	  //{	  
	  pair_pipa = tools::toms748_solve(forceipafunct, min, max, forceipa_rt, forceipa_max, tol, max_iter);

	  //pair_pipa = tools::toms748_solve(forcepipafunct, rt, max, tol, max_iter);
	
	 // std::cerr << "\n" << rt << "\t" << max << "\t" << 
	 //}
	}
	catch(const std::exception& exc) {
	  std::cerr << '\n' << exc.what() << '\n';
	  std::terminate();
	}
	if(max_iter > 990){
	  std::cerr << "\n ERROR max iter " << max_iter << "\n\n";
	  std::terminate();
	}
	double rbarrier = 0.5*(pair_pipa.first+pair_pipa.second);
	if(rbarrier>=0){
	  double phimax = potprefactor * pmpcfunct(rbarrier);
	  double eta = exp(-phimax/(Kboltz*temperature))
	               *(1.0+(phimax-phi_rt)/(Kboltz*temperature));	  
	  return eta;
	}
	else{
	  std::cerr << "\n ERROR Negative rbarrier " << rbarrier << '\n';
	  std::terminate();
	}
      }
    }
    
    return 1.0;
  }
  
}

#endif//EINT_H
