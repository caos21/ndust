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

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <utility>

#include <omp.h>

#include <boost/array.hpp>
#include <boost/math/tools/roots.hpp>
// #include <boost/numeric/odeint.hpp>

#include <boost/math/special_functions/factorials.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_sparse.hpp>
#include <boost/numeric/bindings/lapack/gesv.hpp>
#include <boost/numeric/bindings/traits/ublas_vector2.hpp>

// #include <boost/numeric/bindings/umfpack/umfpack.hpp>

#include "constants.h"

namespace tools = boost::math::tools;
namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;
// namespace umf = boost::numeric::bindings::umfpack;
namespace bmath = boost::math;

namespace eint {
  typedef ublas::matrix<double, ublas::column_major> ubmatrix;
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

  // compute multipolar coefficients
  inline
  void compute_Acoefficients(ubmatrix &coefficients, ubvector &independent,
                            double &A10,
                            const double r1, const double q1,
                            const double r2, const double q2,
                            const double rn,
                            const double eps,
                            const unsigned int j1max=5,
                            const unsigned int j2max=5) {

    A10 = Kcoul*q1;
    for(unsigned int j1=1; j1<j1max; ++j1){
      double prefac1 = (eps-1.0)*j1/((eps+1.0)*j1 + 1.0);
      // term due to Q2, j1!=0, always negative
      double f2 = -prefac1*pow(r1, 2*j1+1)*Kcoul*q2/pow(rn, j1+1);
      // this term is independent of coefficients A
      independent[j1-1] = f2;
      for(unsigned int j3=0; j3<j1max; ++j3){
        for(unsigned int j2=1; j2<j2max; ++j2){
          double prefac2 = (eps-1.0)*j2/((eps+1.0)*j2 + 1.0);
          double denom = (bmath::factorial<double>(j1)
                          *bmath::factorial<double>(j2)
                          *bmath::factorial<double>(j2)
                          *bmath::factorial<double>(j3))
                          *pow(rn, j1+2*j2+j3+2);
          double numer = prefac1 * prefac2 * bmath::factorial<double>(j1+j2)
                      * bmath::factorial<double>(j2+j3)
                      * pow(r1, 2*j1+1) * pow(r2, 2*j2+1);

          if(j3!=0) {
            coefficients(j1-1, j3-1) -= numer / denom;
          }
          else {
            independent[j1-1] += A10*numer / denom;
          }
        }
        if(j3==j1) {
          coefficients(j1-1, j3-1) += 1.0;
        }
      }
    }
  }

  //! Bichoutskaia potential.
  /*!
  * 
  *As defined in equation 9 of 
    Lindgren, E. B., Chan, H.-K., Stace, A. J. & Besley, E.
    Progress in the theory of electrostatic interactions between charged particles.
    Phys. Chem. Chem. Phys. 18, 5883–5895 (2016) (and A7 of Bichoutskaia)
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

    unsigned int size = acoefficients.size();

    double invK = 1.0/Kcoul;

    double pot_coul = Kcoul*q1*q2/h;// Coulomb term

    double pot_2 = 0.0;
    for(unsigned int m=1; m<size; ++m) {// m:1 -> lns
      for(unsigned int l=0; l<size; ++l) {
        double numer = (acoefficients[l]*(eps-1.0)*m
                      *bmath::factorial<double>(l+m)*pow(r2, 2*m+1));
        double denom = ((eps+1.0)*m+1)*(bmath::factorial<double>(l)
                      *bmath::factorial<double>(m))*pow(h, 2*m+l+2);
        pot_2 += numer / denom;
      }
    }
    pot_2 *= -q1;

    double pot_3 = 0.0;
    for(unsigned int l=1; l<size; ++l) {
      double numer = (acoefficients[l]*acoefficients[l])*((eps+1.0)*(l)+1);
      double denom = (eps-1.0)*l*pow(r1, 2*l+1);
      pot_3 += numer / denom;
    }
    pot_3 *= -invK;

    // WARNING 0.5 double count?
    return pot_coul + 0.5*pot_2 + 0.5*pot_3;
  }

  //! Bichoutskaia force.
  /*!
  * 
  *As defined in equation 9 of 
    Lindgren, E. B., Chan, H.-K., Stace, A. J. & Besley, E.
    Progress in the theory of electrostatic interactions between charged particles.
    Phys. Chem. Chem. Phys. 18, 5883–5895 (2016) (and A7 of Bichoutskaia)
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
  double force_bichoutskaia(const double r1, const double q1,
                            const double r2, const double q2,
                            const double h,
                            const ubvector& acoefficients,
                            const double eps) {
    unsigned int size = acoefficients.size();

    double invK = 1.0/Kcoul;

    double force_coul = Kcoul*q1*q2/(h*h);// Coulomb term

    double force_2 = 0.0;
    for(unsigned int m=1; m<size; ++m) {// m:1 -> lns
      for(unsigned int l=0; l<size; ++l) {
        double numer = (acoefficients[l]*(eps-1.0)*m*(m+1)
                      *bmath::factorial<double>(l+m)*pow(r2, 2*m+1));
        double denom = ((eps+1.0)*m+1)*(bmath::factorial<double>(l)
                      *bmath::factorial<double>(m))*pow(h, 2*m+l+3);
        force_2 += numer / denom;
      }
    }
    force_2 *= -q1;

    double force_3 = 0.0;
    for(unsigned int l=1; l<size-1; ++l) {
      double numer = (acoefficients[l]*acoefficients[l+1])*((eps+1.0)*(l+1)+1);
      double denom = (eps-1.0)*pow(r1, 2*l+3);
      force_3 += numer / denom;
    }
    force_3 *= -invK;

    return force_coul + force_2 + force_3;
  }

  inline
  double relative_error(const double trueval, const double expval){
    return abs((expval-trueval)/trueval);
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

  // compute ipa enhancement factor for a pair of particles
  inline
  double efactor_ipa(double r1, double r2,
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
    
    potential_ipa_funct pipafunct(r21, q21, eps);
    double pipa = pipafunct(rt);

    std::pair<double, double> pair_pipa;
    bool failed = false;
    try {
      pair_pipa = tools::toms748_solve(pipafunct, min, max, tol, max_iter);
    }
    catch(const std::exception& e) {
      failed = true;
      pair_pipa.first = 0.0;
    }
    bool attcontact = (pipa<0.0? true: false);
    bool fullatt = attcontact && failed;
    bool withphimax = !failed || attcontact;

    double potprefactor = q1*q1*eCharge*eCharge/r1;
    double eta = 0.0;
    double phimin = potprefactor*pipa;

    if(withphimax){
      double phimax = potprefactor*pair_pipa.first;
      eta = exp(-phimax/(Kboltz*temperature))
          *(1.0+(phimax-phimin)/(Kboltz*temperature));
    }
    if(fullatt){
      eta = 1.0 - phimin /(Kboltz*temperature);
    }
    if(!attcontact){
      eta = exp(-phimin/(Kboltz*temperature));
    }

  //   std::cout << std::endl << r21 << "\t" << q21 << "\t" << potprefactor*pair_pmpc.first
  //             << "\t" << phimin << "\t" << attcontact << "\t"
  //             << fullatt << "\t" << eta;
    return eta;
  }


  // multipolar coefficients potential functor
  struct potential_mpc_funct
  {
    potential_mpc_funct(double r21_,
                        double q21_,
                        double eps_,
                        unsigned int j1max_,
                        unsigned int j2max_):
                        r21(r21_), q21(q21_), eps(eps_),
                        j1max(j1max_), j2max(j2max_){

      coefficients = ublas::zero_matrix<double>(j1max-1, j1max-1);
      independent = ublas::zero_vector<double>(j1max-1);
      acoefficients = ublas::zero_vector<double>(j1max);
    }

    double operator()(double const& rt) {
      double A10=0.0;

      // compute coefficients
      compute_Acoefficients(coefficients, independent, A10, 1.0, 1.0,
                            r21, q21, rt, eps, j1max, j2max);
      // solve linear system
      lapack::gesv(coefficients, independent);
      
      // fill acoefficients
      acoefficients[0] = A10;
      for(unsigned int l=0; l<independent.size(); ++l){
        acoefficients[l+1] = independent[l];
      }

      // compute Bichoutskaia potential
      return potential_bichoutskaia(1.0, 1.0, r21, q21, rt,
                                    acoefficients, eps);

    }

    double r21;
    double q21;
    double eps;
    unsigned int j1max;
    unsigned int j2max;
    ubmatrix coefficients;
    ubvector independent;
    ubvector acoefficients;
  };

  // compute eta factor for a pair of particles
  inline
  double efactor_mpc(double r1, double r2,
                    double q1, double q2,
                    double eps, double temperature,
                    unsigned int j1max, unsigned int j2max) {
    
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
    
    potential_mpc_funct pmpcfunct(r21, q21, eps, 
                                  j1max, j2max);
    double pmpc = pmpcfunct(rt);

    std::pair<double, double> pair_pmpc;
    bool failed = false;
    try {
      pair_pmpc = tools::toms748_solve(pmpcfunct, min, max, tol, max_iter);
    }
    catch(const std::exception& e) {
      failed = true;
      pair_pmpc.first = 0.0;
    }
    bool attcontact = (pmpc<0.0? true: false);
    bool fullatt = attcontact && failed;
    bool withphimax = !failed || attcontact;

    double potprefactor = q1*q1*eCharge*eCharge/r1;
    double eta = 0.0;
    double phimin = potprefactor*pmpc;

    if(withphimax){
      double phimax = potprefactor*pair_pmpc.first;
      eta = exp(-phimax/(Kboltz*temperature))
          *(1.0+(phimax-phimin)/(Kboltz*temperature));
    }
    if(fullatt){
      eta = 1.0 - phimin /(Kboltz*temperature);
    }
    if(!attcontact){
      eta = exp(-phimin/(Kboltz*temperature));
    }

  //   std::cout << std::endl << r21 << "\t" << q21 << "\t" << potprefactor*pair_pmpc.first
  //             << "\t" << phimin << "\t" << attcontact << "\t"
  //             << fullatt << "\t" << eta;
    return eta;
  }
}

#endif//EINT_H
