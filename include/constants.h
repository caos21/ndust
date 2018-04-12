/*
 * Copyright 2016 Benjamin Santos <caos21@gmail.com>
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

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>


// math constants, pi
#include <boost/math/constants/constants.hpp>

#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_const_num.h>

//Relative error for element in profile
// #define PROF_REL_ERR 0.05//0.0005 optimal
#define EPSILONZETA 1.0e-250
#define EPSILON 1.0e-150
// #define EPSILONETA 1.0e-50
#define EPSILONETA 1.0e-18
// #define EPSILONDEATH 1.0e-50
#define EPSILONDEATH 1.0e-18//1.0e-30
//
#define ZETA1_12 0.0833333333333
#define ZETA1_720 0.00138888888889

// factorials up to 90!
const double facts[]= {
  1.000000000000000000e+00,
  1.000000000000000000e+00,
  2.000000000000000000e+00,
  6.000000000000000000e+00,
  2.400000000000000000e+01,
  1.200000000000000000e+02,
  7.200000000000000000e+02,
  5.040000000000000000e+03,
  4.032000000000000000e+04,
  3.628800000000000000e+05,
  3.628800000000000000e+06,
  3.991680000000000000e+07,
  4.790016000000000000e+08,
  6.227020800000000000e+09,
  8.717829120000000000e+10,
  1.307674368000000000e+12,
  2.092278988800000000e+13,
  3.556874280960000000e+14,
  6.402373705728000000e+15,
  1.216451004088320000e+17,
  2.432902008176640000e+18,
  5.109094217170944000e+19,
  1.124000727777607680e+21,
  2.585201673888497821e+22,
  6.204484017332394100e+23,
  1.551121004333098391e+25,
  4.032914611266057190e+26,
  1.088886945041835194e+28,
  3.048883446117138016e+29,
  8.841761993739700773e+30,
  2.652528598121911042e+32,
  8.222838654177923583e+33,
  2.631308369336935547e+35,
  8.683317618811887119e+36,
  2.952327990396041951e+38,
  1.033314796638614664e+40,
  3.719933267899011775e+41,
  1.376375309122634558e+43,
  5.230226174666011171e+44,
  2.039788208119744159e+46,
  8.159152832478980083e+47,
  3.345252661316381834e+49,
  1.405006117752880287e+51,
  6.041526306337383407e+52,
  2.658271574788448529e+54,
  1.196222208654801886e+56,
  5.502622159812089850e+57,
  2.586232415111681777e+59,
  1.241391559253607253e+61,
  6.082818640342675225e+62,
  3.041409320171338142e+64,
  1.551118753287382377e+66,
  8.065817517094387685e+67,
  4.274883284060025485e+69,
  2.308436973392413792e+71,
  1.269640335365827802e+73,
  7.109985878048635815e+74,
  4.052691950487721410e+76,
  2.350561331282878495e+78,
  1.386831185456898386e+80,
  8.320987112741389895e+81,
  5.075802138772248358e+83,
  3.146997326038793939e+85,
  1.982608315404440306e+87,
  1.268869321858841513e+89,
  8.247650592082471517e+90,
  5.443449390774430694e+92,
  3.647111091818868322e+94,
  2.480035542436831022e+96,
  1.711224524281412974e+98,
  1.197857166996989027e+100,
  8.504785885678624248e+101,
  6.123445837688610055e+103,
  4.470115461512684895e+105,
  3.307885441519386893e+107,
  2.480914081139539975e+109,
  1.885494701666050381e+111,
  1.451830920282858721e+113,
  1.132428117820629898e+115,
  8.946182130782977113e+116,
  7.156945704626380571e+118,
  5.797126020747367841e+120,
  4.753643337012842020e+122,
  3.945523969720658788e+124,
  3.314240134565353194e+126,
  2.817104114380551285e+128,
  2.422709538367273905e+130,
  2.107757298379527854e+132,
  1.854826422573984361e+134,
  1.650795516090846024e+136
};

// pi from boost
const double pi = boost::math::double_constants::pi;

// TODO use boost

const double Kboltz = GSL_CONST_MKSA_BOLTZMANN;//  kg m^2 / K s^2 
// const double Pressure = 13.332237;
// const double Temperature = 300.0;
const double EpsilonZero = GSL_CONST_MKSA_VACUUM_PERMITTIVITY;
const double Kcoul = 1.0/(4.0*M_PI*EpsilonZero);
const double eCharge = GSL_CONST_MKSA_ELECTRON_CHARGE;
const double Kboltz_eV = Kboltz / eCharge; // eV / K
const double eMass = GSL_CONST_MKSA_MASS_ELECTRON;
// const double NAvogadro = GSL_CONST_NUM_AVOGADRO; // 1/mol
// const double ArMMass = 39.948; // Argon Molecular Mass (NIST) in g/mol
// const double ArMass = (ArMMass/NAvogadro)/1000.0; // Argon Mass in Kg
// const double Acc = 0.9; // Accommodation factor alpha_m
// const double SiDensity = 2329.0; //Kg/m3 https://en.wikipedia.org/wiki/Silicon
// const double EpsilonSi = 11.68; // 11.7 in
// const double AffinitySi = 4.05;//eV
// http://www.pveducation.org/pvcdrom/materials/general-properties-of-silicon

// const double critical_diameter = 1.8E-9; // Multiple charges in particle

// const double maximum_charge = 1.0E4; // Maximum charge 
#endif
