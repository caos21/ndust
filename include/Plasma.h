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

#ifndef PLASMA_H
#define PLASMA_H

const boost::uintmax_t EPS_MAXITER = 10000;
const boost::math::tools::eps_tolerance<double> EPS_TOL(30);

inline
double temperature_from_energy(const double energy_ev) {
	return 2.0*energy_ev/(3.0*Kboltz_eV);
}

inline
double arrhenius_rate(const double energy_ev, const double activation_energy,
					  const double const_A, const double const_beta) {
		double temperature = temperature_from_energy(energy_ev);
		return const_A*pow(temperature, const_beta)
				*exp(-3.0*activation_energy/(2.0*energy_ev));
}

//
// converts from cm^3 / mol.s -> m^3/s
// # A in cm3 /s mol  x 1e-6 m3/cm3 / Na(mol^-1) -> [A] in m3/s
//
inline
double aconst_convert(const double const_A) {
	return const_A * 1.e-6/NAvogadro;
}

//
// converts from cal/mol to eV
// # Ea in cal/mol * 4.184 J/cal * (1/1.6x10-19)eV/J / Na(mol^-1) -> Ea in eV
//
inline
double aenergy_convert(const double activation_energy) {
    return activation_energy * 4.184 / (eCharge*NAvogadro);
}

inline
std::pair<double, double> arrhenius_values(const double const_A,
										   const double activation_energy) {
	return std::make_pair(aconst_convert(const_A),
						  aenergy_convert(activation_energy));
}

inline
std::pair<double, double> arrhenius_values(std::pair<double, double> aconv) {
	return arrhenius_values(aconv.first, aconv.second);
}

/**
 * @brief Represents rates
 * 
 */
class Rate {
public:

	Rate(): rate_name(""),
		 activation_energy(0.0),
		 const_A(0.0),
		 const_beta(0.0) {

	}

	Rate(std::string rate_name, double activation_energy,
	     double const_A, double const_beta=0.0):
		 rate_name(rate_name),
		 activation_energy(activation_energy),
		 const_A(const_A),
		 const_beta(const_beta) {
	}

	double operator() (const double energy_ev) const {
		return arrhenius_rate(energy_ev, activation_energy,
							  const_A, const_beta);
	}

	void print() const {
		std::cerr << std::endl << "----------------------------------------------------"
				  << std::endl;
		std::cerr << "Rate: " << rate_name << std::endl;
		std::cerr << "E_A: "<< activation_energy << std::endl;
		std::cerr << "A: " << const_A << std::endl;
		std::cerr << "beta: " << const_beta << std::endl;
	}
	
	std::string rate_name;
	double activation_energy;
	double const_A;
	double const_beta; 
};

typedef std::map<std::string, std::shared_ptr<Rate> > ratemap;

class Plasma {
public:
	Plasma() {}

	Plasma(ratemap& rates_map, double ni=1.1e16, double meta_on=1.0):
            rates_map(rates_map), meta_on(meta_on) {
		
	    ion_dens = std::make_pair(ni, ni);
		e_dens = std::make_pair(ni, ni);
		meta_dens = std::make_pair(ni, ni);
	}

	std::pair<double, double> operator()(double const& energy) { 
    	return std::make_pair(full_ionbalance(energy),
							  diff_ionbalance(energy));
  	}

	// WARNING functor for toms and bracket and solve
	// double operator()(double const& energy) {
    // 	return full_ionbalance(energy);
	// }
    // double operator()(double const& energy) {
    //     if(meta_on>0.0) {
    //         return full_ionbalance(energy);
    //     }
    //     else {
    //         return ionbalance(energy);
    //     }
    // }

	double ambipolar_diffusion(double energy) {
		double ne = e_dens.first;
		double ni = ion_dens.first;
        double ng = gas_dens;
		// return (ion_diffusion*((1.+2.*energy/(3.*Kboltz_eV * temperature))
		// 		/(1.+(ni/ne)*(ion_mobility/e_mobility))));
        double alpha = ng/ne;
        double gamma = 2.*energy/(3.*Kboltz_eV * temperature);
        return ion_diffusion * 1.0 + gamma + 2.0*alpha*gamma / (1.0 + alpha*gamma);
	}

    double ionbalance(double energy) {
		Rate &ki = *rates_map["ionization"];

		double ng = gas_dens;
		double ne = e_dens.first;
		double ni = ion_dens.first;

		return ki(energy)*ng*ne
			  -ni*ambipolar_diffusion(energy)/pow(length/pi, 2)
              + 2.8e-14*ni*ne
              -ion_loss.first;
	}

	double full_ionbalance(double energy) {
		Rate &ki = *rates_map["ionization"];
		Rate &ksi = *rates_map["stepwise_ionization"];
		Rate &kmp = *rates_map["metastable_pooling"];

		double ng = gas_dens;
		double ne = e_dens.first;
		double ni = ion_dens.first;
		double nm = meta_dens.first;

		return ki(energy)*ng*ne
			  +ksi(energy)*nm*ne
			  +kmp(energy)*nm*nm
			  -ni*ambipolar_diffusion(energy)/pow(length/pi, 2)
              -ion_loss.first;
	}

	double diff_ionbalance(double energy) {
		Rate &ki = *rates_map["ionization"];
		Rate &ksi = *rates_map["stepwise_ionization"];
		
		double ng = gas_dens;
		double ne = e_dens.first;
		double ni = ion_dens.first;
		double nm = meta_dens.first;

		double ion_difffactor = (ki.const_beta/energy
								 +(3.0*ki.activation_energy/
								  (2.0*energy*energy)));
		
		double ksi_difffactor = (3.0*ksi.activation_energy/
								 (2.0*energy*energy));

		return ki(energy)*ng*ne*ion_difffactor
			  +ksi(energy)*nm*ne*ksi_difffactor
			  -ni*2.0*ion_diffusion/(3.0*Kboltz_eV * temperature*pow(length/pi, 2)*(1.0+(ion_mobility/e_mobility)));
	}

	int meta_evolution(double dt) {
		Rate &kex = *rates_map["excitation"];
		Rate &ksi = *rates_map["stepwise_ionization"];
		Rate &ksc = *rates_map["superelastic_colision"];
		Rate &kmp = *rates_map["metastable_pooling"];

		Rate &kr = *rates_map["quenching_toresonant"];
		Rate &k2q = *rates_map["twobody_quenching"];
		Rate &k3q = *rates_map["threebody_quenching"];

		double energy = mean_energy.second;
		double ng = gas_dens;
		double ne = e_dens.first;
		double nm = meta_dens.first;

		meta_dens.second = nm + dt
		  * (kex(energy)*ng*ne - ksi(energy)*nm*ne - ksc(energy)*nm*ne
		    - kr(energy)*nm*ne - 2.0*kmp(energy)*nm*nm - k2q(energy)*ng*nm
		 	- k3q(energy)*ng*ng*nm - nm*meta_diffusion*pow(pi/length, 2));

		return 0;
	}

    int electron_evolution(double dt, double nano_qdens=0.0) {
        // Rate &ki = *rates_map["ionization"];
		// Rate &ksi = *rates_map["stepwise_ionization"];
		// Rate &kmp = *rates_map["metastable_pooling"];

        // double energy = mean_energy.second;
		// double ng = gas_dens;
		// double ne = e_dens.first;
		// double nm = meta_dens.first;

        // e_dens.second = ne + dt
        //     * (ki(energy)*ng*ne + ksi(energy)*nm*ne + kmp(energy)*nm*nm
        //        -e_loss);
    
        e_dens.second = ion_dens.first - nano_qdens;
        return 0;
    }
	
	std::pair<double, double> mean_energy;
	std::pair<double, double> ion_dens;
	std::pair<double, double> e_dens;
	std::pair<double, double> meta_dens;

	std::pair<double, double> ion_loss;

    double meta_on;

private:
	double length = 3e-2;// m
	double temperature = 300.0;
	double ion_temperature = 300.0;
	double KbTg = Kboltz * temperature;
	double pressure_torr = 0.1;
	double pressure = 133.32237 * pressure_torr;

	double gas_dens = pressure / KbTg;

	double ion_diffusion = 4.0e-3 / pressure_torr;// m2/s	
	double ion_mobility = 0.1444 / pressure_torr;// m2/V.s

	double electron_diffusion = 120.0 / pressure_torr;// m2/s
	double e_mobility = 30.0 / pressure_torr;// m2/V.s

	double meta_diffusion = 2.42e20/gas_dens;
	ratemap rates_map;
};


#endif// PLASMA_H