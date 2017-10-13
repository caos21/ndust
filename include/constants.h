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

#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_const_num.h>

//Relative error for element in profile
// #define PROF_REL_ERR 0.05//0.0005 optimal
#define EPSILONZETA 1.0e-250
#define EPSILON 1.0e-150
// #define EPSILONETA 1.0e-50
#define EPSILONETA 1.0e-50
// #define EPSILONDEATH 1.0e-50
#define EPSILONDEATH 1.0e-50//1.0e-30
//
#define ZETA1_12 0.0833333333333
#define ZETA1_720 0.00138888888889

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
