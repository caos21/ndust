/*
 * Copyright 2019 <Benjamin Santos> <caos21@gmail.com>
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

#ifndef PLASMACHEM_H
#define PLASMACHEM_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <memory>
#include <algorithm>
#include <functional>
#include <cmath>
#include <tuple>
#include <array>
#include <map>

// hdf5 c++ bindings
#include <H5Cpp.h>

#include "constants.h"

extern "C" {
#include "common.h"
#include "lsoda.h"
#include "lsoda_internal.h"
#include "blas.h" 
}


/**
 * @brief Represents rates
 * 
 */
class RateSpec {
public:

	RateSpec(): rate_function(""),
		avar(0.0),
		bvar(0.0),
		cvar(0.0) {
		set_params();
	}

	RateSpec(std::string rate_function,
		double avar,
		double bvar = 0.0,
		double cvar = 0.0):
		rate_function(rate_function),
		avar(avar),
		bvar(bvar),
		cvar(cvar) {
		set_params();
	}

	void set_params() {
		params[0] = avar;
		params[1] = bvar;
		params[2] = cvar;
	}

	double operator() (const double energy_ev) const {
		if (rate_function == "arrhenius") {
			return avar*pow(energy_ev, cvar)*exp(-bvar/energy_ev);
			//arrhenius_rate(energy_ev, avar, bvar, cvar);
		}
		if (rate_function == "constant") {
			return avar;
		}
		if (rate_function == "a1expb") {
			return avar * (1.0-exp(-bvar*energy_ev));
		}
		return 0.0;
	}

	void print() const;
	
	std::string rate_function;
	double avar;
	double bvar;
	double cvar;
	double params[3];
};

typedef std::map<std::string, std::shared_ptr<RateSpec> > Rates;

class PlasmaChem {
public:
	PlasmaChem() {}

	PlasmaChem(Rates& rates_map):
              rates_map(rates_map) {
	}

	void init_rates();

	void init_parameters(double length = 4e-2,
						 double radius = 6e-2,
						 double temperature = 300.0,
						 double ion_temperature = 300.0,
						 double pressure_torr = 0.1,
					     double arsih4_ratio = 29.0/30.0,
						 double armass =   6.6335209e-26,
						 double sih4mass = 5.3150534e-26,
						 double power = 1.0);

	// WARNING useless
	void update_density(std::vector<double> pdens, std::vector<double> ndens);

	void update_nSiH4(double density) {
		nSiH4 = density;
	}

	int solve(const double ti,
			  const double tf,
			  const double dt,
			  const std::vector<double> &density_sourcedrain,
			  const double energy_sourcedrain,
			  const std::vector<double> &olddensity,
			  std::vector<double> &newdensity,
			  double &min_dtq);

	double thermal_velocity(double mass);

	double diffusion_neutrals(double mass, double lamda=3.5*1e-3);

 	double center2edge_neutrals(double mass);

	double flux_neutrals(double mass);
    
	double bohm_velocity(double mass);
    
	double center2edge_ions(double mass, double lambda);

	double flux_ions(double mass, double lambda);

	double ion_velocity(double mass);
    
	friend int plasma_system(double t, double *n, double *dndt, void *data);
	
	std::vector<double> ndens;
	std::vector<double> pdens;
	std::vector<double> density_sourcedrain;
	
	double nano_qdens;
	double nano_qdens_rate;
	//std::pair<double, double> ion_dens;

	Rates rates_map;
	//std::shared_ptr<Rates> rates_map = std::make_shared<Rates>();

private:
	double length;
	double radius;
	double reactor_volume;
	double reactor_area;
	double ratio_AV;
	double temperature;
	double ion_temperature;
	double KbTg;
	double pressure_torr;
	double pressure;
	double arsih4_ratio ;
	double armass;
	double sih4mass;
	double power;
	double gas_dens;
	double nAr;
	double nSiH4;
	double vth_Ar;
	double vth_SiH4;
	double flux_SiH3;
	double flux_SiH2;
	double flux_Ar;
	double flux_Arp;
	double flux_SiH3p;
	double lambdai;
	double vsheath;

	// double ion_diffusion = 4.0e-3 / pressure_torr;// m2/s	
	// double ion_mobility = 0.1444 / pressure_torr;// m2/V.s

	// double electron_diffusion = 120.0 / pressure_torr;// m2/s
	// double e_mobility = 30.0 / pressure_torr;// m2/V.sx

	// double meta_diffusion = 2.42e20/gas_dens;
	
};

inline
int plasma_system(double t, double *n, double *dndt, void *data) {

  	PlasmaChem* plasmachem = static_cast< PlasmaChem* >(data);

	// std::shared_ptr<RateSpec> ki = plasmachem->rates_map["R2:ki"];	
	// n[0] = (*ki)(t);
	
	double ne     = n[0];
    double nArp   = n[1];
    double nArm   = n[2];
    double nSiH3p = n[3];
    double nSiH3  = n[4];
    double nSiH2  = n[5];
    double neps   = n[6];

	double energy = neps/ne;
    
	double kel = plasmachem->rates_map["R1:kel"]->operator()(energy);
	double ki  = plasmachem->rates_map["R2:ki"]->operator()(energy);
	double kex = plasmachem->rates_map["R3:kex"]->operator()(energy);
	double kelsih4 = plasmachem->rates_map["R4:kelsih4"]->operator()(energy);
	double kdisih4 = plasmachem->rates_map["R5:kdisih4"]->operator()(energy);
	double kdsih3 = plasmachem->rates_map["R6:kdsih3"]->operator()(energy);
	double kdsih2 = plasmachem->rates_map["R7:kdsih2"]->operator()(energy);
	double kisih3 = plasmachem->rates_map["R8:kisih3"]->operator()(energy);
	double kv13 = plasmachem->rates_map["R9:kv13"]->operator()(energy);
	double kv24 = plasmachem->rates_map["R10:kv24"]->operator()(energy);
	double k12 = plasmachem->rates_map["R12:k12"]->operator()(energy);
	double k13 = plasmachem->rates_map["R13:k13"]->operator()(energy);
	double k14 = plasmachem->rates_map["R14:k14"]->operator()(energy);
	double k15 = plasmachem->rates_map["R15:k15"]->operator()(energy);	

	double eki = plasmachem->rates_map["R2:ki"]->bvar;
	double ekex = plasmachem->rates_map["R3:kex"]->bvar;
	double ekdisih4 = plasmachem->rates_map["R5:kdisih4"]->bvar;
	double ekdsih3 = plasmachem->rates_map["R6:kdsih3"]->bvar;
	double ekdsih2 = plasmachem->rates_map["R7:kdsih2"]->bvar;
	double ekisih3 = plasmachem->rates_map["R8:kisih3"]->bvar;
	double ekv13 = plasmachem->rates_map["R9:kv13"]->bvar;
	double ekv24 = plasmachem->rates_map["R10:kv24"]->bvar;

	double nAr = plasmachem->nAr;
	double nSiH4 = plasmachem->nSiH4;
	double flux_Arp = plasmachem->flux_Arp;
	double flux_Ar = plasmachem->flux_Ar;
	double flux_SiH3p = plasmachem->flux_SiH3p;
	double flux_SiH3 = plasmachem->flux_SiH3;
	double flux_SiH2 = plasmachem->flux_SiH2;
	double ratio_AV = plasmachem->ratio_AV;
	double volume = plasmachem->reactor_volume;

	std::vector<double> sourcedrain = plasmachem->density_sourcedrain;

	//nArp = ne - nSiH3p + plasmachem->nano_qdens;

	// WORKING
    double dne = +ki*nAr*ne + kdisih4*ne*nSiH4 + kisih3*ne*nSiH3
			 //- flux_Arp*ratio_AV*nArp - flux_SiH3p*ratio_AV*nSiH3p
			 /*- sourcedrain[1]*volume*/ - sourcedrain[0]*volume;// + sourcedrain[0] <- workin with this term;
    //double dnArp = dne - dSiH3p;
	//

	// double dne = +ki*nAr*ne + kdisih4*ne*nSiH4 + kisih3*ne*nSiH3
	// 		 - flux_Arp*ratio_AV*nArp - flux_SiH3p*ratio_AV*nSiH3p
	// 		 - sourcedrain[1]*volume - sourcedrain[0]*volume;// + sourcedrain[0] <- workin with this term;

	double dnArp  = +ki*nAr*ne - flux_Arp*ratio_AV*nArp
	 			  - sourcedrain[1]*volume;// - sourcedrain[0]*volume ;// + sourcedrain[0];

	// double dne = +ki*nAr*ne + kdisih4*ne*nSiH4 + kisih3*ne*nSiH3
	// 		 - flux_Arp*ratio_AV*nArp - flux_SiH3p*ratio_AV*nSiH3p
	// 		 - sourcedrain[1];
    // double dnArp  = +ki*nAr*ne - flux_Arp*ratio_AV*nArp
	// 			  - sourcedrain[0];
	//dnArp =0.0;
    double dnArm  = +kex*nAr*ne -k12*nArm*nSiH4 -k13*nArm*nSiH4 -k14*nArm*nSiH3 -k15*nArm*nSiH2 -flux_Ar*ratio_AV*nArm;
    //double #dSiH4 = -kdisih4*ne*nSiH4 -kdsih3*ne*nSiH4 -kdsih2*ne*nSiH4 -kv13*ne*nSiH4 -kv24*ne*nSiH4 
    double dSiH3p = +kdisih4*ne*nSiH4 + kisih3*ne*nSiH3 - flux_SiH3p*ratio_AV*nSiH3p;
    double dSiH3  = +kdsih3*ne*nSiH4  -kisih3*ne*nSiH3 +k12*nArm*nSiH4 -k14*nArm*nSiH3
				  - flux_SiH3*ratio_AV*nSiH3;
    double dSiH2  = +kdsih2*ne*nSiH4 +k13*nArm*nSiH4 +k14*nArm*nSiH3 -k15*nArm*nSiH2
	              - flux_SiH2*ratio_AV*nSiH2;
    
    //double #dne = dnArp + dSiH3p
	//double dnArp = dne - dSiH3p;
    
	double power = plasmachem->power;
	double reactor_volume = plasmachem->reactor_volume;
	double vsheath = plasmachem->vsheath;
	double emass = eMass;
	double e = eCharge;
	double armass = plasmachem->armass;
	double sih4mass = plasmachem->sih4mass;

    double dneps = (power/reactor_volume
             - eki*ki*nAr*ne
             - ekex*kex*nAr*ne
             - (5./3.)*plasmachem->bohm_velocity(armass)*ratio_AV*neps
             - e*vsheath*plasmachem->bohm_velocity(armass)*ratio_AV*ne
             - 3.0*(emass/armass)*kel*neps*nAr
             - 3.0*(emass/sih4mass)*kelsih4*neps*nSiH4
             - ekisih3*kisih3*ne*nSiH3
             - ekdisih4*kdisih4*ne*nSiH4
             - ekdsih3*kdsih3*ne*nSiH4
             - ekdsih2*kdsih2*ne*nSiH4
             - ekv13*kv13*ne*nSiH4
             - ekv24*kv24*ne*nSiH4
			 - sourcedrain[6]*volume);
    
    
	dndt[0] = dne==dne ? dne : 0.0;
	dndt[1] = dnArp==dnArp? dnArp : 0.0;
	dndt[2] = dnArm==dnArm? dnArm : 0.0;
	dndt[3] = dSiH3p==dSiH3p? dSiH3p : 0.0;
	dndt[4] = dSiH3==dSiH3? dSiH3 : 0.0;
	dndt[5] = dSiH2==dSiH2? dSiH2 : 0.0;
	dndt[6] = dneps==dneps? dneps : 0.0;

	return 0;
}


// WARNING lorenz system test
inline
int lorenz_system(double t, double *n, double *dndt, void *data) {

  	PlasmaChem* plasmachem = static_cast< PlasmaChem* >(data);

	double x = n[0];
	double y = n[1];
	double z = n[2];

	double rho = plasmachem->rates_map["lorenz"]->avar;
	double sigma = plasmachem->rates_map["lorenz"]->bvar;
	double beta = plasmachem->rates_map["lorenz"]->cvar;
	
	dndt[0] = sigma*(y-x);
  	dndt[1] = rho*x - y - x*z;
  	dndt[2] = x*y - beta*z;
	
	return 0;
}
#endif// PLASMACHEM_H
