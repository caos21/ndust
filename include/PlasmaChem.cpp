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

#include "PlasmaChem.h"

/**
 * @brief Represents rates
 * 
 */

void RateSpec::print() const {
    std::cerr << std::endl << "----------------------------------------------------"
                << std::endl;
    std::cerr << "Rate function: " << rate_function << std::endl;
    std::cerr << "E_A: "<< avar << std::endl;
    std::cerr << "A: " << bvar << std::endl;
    std::cerr << "beta: " << cvar << std::endl;
}

void PlasmaChem::init_rates() {
    rates_map["R1:kel"] = std::make_shared<RateSpec>(RateSpec("a1expb",
                                                        1.8560282921521678e-13,
                                                        0.1633129254054345,
                                                        0.0));
    rates_map["R2:ki"] = std::make_shared<RateSpec>(RateSpec("arrhenius",
                                                                7.41485579e-14,
                                                                15.8,
                                                                5.00954185e-3));
    rates_map["R3:kex"] = std::make_shared<RateSpec>(RateSpec("arrhenius",
                                                    3.364506180166679e-14,
                                                    11.5,
                                                    3.7596917572476254e-09));
    // rates_map["stepwise_ionization"]  = std::make_shared<RateSpec>(
    //                         RateSpec("arrhenius", 7.52e16, 124191.0, 0.00));	
    // rates_map["superelastic_colision"] = std::make_shared<RateSpec>(
    //                     RateSpec("arrhenius", 2.60e14, 0.0, 0.74));
    // rates_map["metastable_pooling"] = std::make_shared<RateSpec>(
    //                     RateSpec("arrhenius", 3.70e14, 0.0, 0.0));
    // in m3/s
    // rates_map["quenching_toresonant"] = std::make_shared<RateSpec>(
    //                                         RateSpec("constant", 2.0e-13));
    // rates_map["twobody_quenching"] = std::make_shared<RateSpec>(
    //                                         RateSpec("constant", 3.0e-21));
    // rates_map["threebody_quenching"] = std::make_shared<RateSpec>(
    //                                         RateSpec("constant", 1.1e-43));
    rates_map["R4:kelsih4"] = std::make_shared<RateSpec>(RateSpec("a1expb",
                                                        2.663032076381496e-13,
                                                        0.27729017757606,
                                                        0.0));
    rates_map["R5:kdisih4"] = std::make_shared<RateSpec>(RateSpec("arrhenius",
                                                               3.6864532111053647e-14,
                                                               12.3,
                                                               7.585426699230944e-09));
    rates_map["R6:kdsih3"] = std::make_shared<RateSpec>(RateSpec("arrhenius",
                                                               1.060628946494521e-13,
                                                               8.4,
                                                               0.07797295284736888));
    rates_map["R7:kdsih2"] = std::make_shared<RateSpec>(RateSpec("arrhenius",
                                                               8.999280746645033e-14,
                                                               7.7,
                                                               9.918880635931725e-10));

    rates_map["R8:kisih3"] = std::make_shared<RateSpec>(RateSpec("arrhenius",
                                                               4.6999933119522564e-14,
                                                               8.0,
                                                               0.1625740097963909));
    rates_map["R9:kv13"] = std::make_shared<RateSpec>(RateSpec("arrhenius",
                                                               5.2570315025458126e-15,
                                                               0.27,
                                                               7.675878736753314e-09));
    rates_map["R10:kv24"] = std::make_shared<RateSpec>(RateSpec("arrhenius",
                                                               9.646752429918597e-15,
                                                               0.113,
                                                               1.5565618639018595e-08));
    rates_map["R12:k12"] = std::make_shared<RateSpec>(RateSpec("constant",
                                                                1.4e-16,
                                                                0.0,
                                                                0.0));
    rates_map["R13:k13"] = std::make_shared<RateSpec>(RateSpec("constant",
                                                                2.59e-16,
                                                                0.0,
                                                                0.0));
    rates_map["R14:k14"] = std::make_shared<RateSpec>(RateSpec("constant",
                                                                9.96e-17,
                                                                0.0,
                                                                0.0));
    rates_map["R15:k15"] = std::make_shared<RateSpec>(RateSpec("constant",
                                                                9.96e-17,
                                                                0.0,
                                                                0.0));    
    std::stringstream ss;

    ss << std::setw(24) << "Rate";
    ss << std::setw(24) << "Constants";
    // std::cout << std::endl <<  << ss.str();
    // std::cout << std::endl <<  << "---------------------------------"
    // 							"---------------------------------";  
    for (auto m : rates_map) {
        ss.str("");
        ss << std::setw(24);
        //ss << m.first;
        ss << std::endl << m.first << '\t' << m.second->rate_function;
        for (auto v: m.second->params) {
            ss << std::setw(12);
            ss << v;
        }
        // std::cout << std::endl <<  << ss.str();
        std::cout << ss.str();
    }
    std::cout << std::endl;
    std::cout << std::endl;

    // ratemap rates_map;

    // for (auto m : rates) {
    // 	std::shared_ptr<Rate> r(new Rate(m.first, m.second[1], m.second[0], m.second[2]));
    // 	rates_map[m.first] = r;
    // }

    double epsi = 0.1;
    double epsf = 25.0;
    size_t N = 100;
    double step = (epsf-epsi) / (N-1);

    std::fstream *rate_file = 
                    new std::fstream("init_rates.dat",
                                    std::fstream::out);

    double eps = epsi;
    for (size_t i=0; i<N; ++i) {
        *rate_file << eps;
        for (auto mr : rates_map) {
            RateSpec &r = *mr.second;
            *rate_file << '\t' << r(eps);
        }
        *rate_file << std::endl;
        eps += step;
    }

    rate_file->close();
    delete(rate_file);
}

void PlasmaChem::init_parameters(double length,
                                 double radius,
						         double temperature,
						         double ion_temperature,
						         double pressure_torr,
					             double arsih4_ratio,
                                 double armass,
                                 double sih4mass,
                                 double power
                                 ) {

    this->length = length;
    this->radius = radius;
    this->temperature = temperature;
    this->ion_temperature = ion_temperature;
    this->pressure_torr = pressure_torr;
    this->arsih4_ratio = arsih4_ratio;
    this->armass = armass;
    this->sih4mass = sih4mass;
    this->power = power/eCharge;//eV/s

    nano_qdens = 0.0;
    nano_qdens_rate = 0.0;
    KbTg = Kboltz * temperature;
	pressure = 133.32237 * pressure_torr;
    reactor_volume = length*pi*radius*radius;
    reactor_area = length*2.0*pi*radius;
    ratio_AV = reactor_area / reactor_volume;
	gas_dens = pressure / KbTg;
	nAr = arsih4_ratio * gas_dens;
	nSiH4 = (1.0-arsih4_ratio) * gas_dens;
    vth_Ar = thermal_velocity(armass);
    vth_SiH4 = thermal_velocity(sih4mass);

    flux_SiH3 = flux_neutrals(sih4mass);
    flux_SiH2 = flux_neutrals(sih4mass);
    flux_Ar = flux_neutrals(armass);

    //From Lieberman pag 80 (117)
    lambdai = 1. / (330 * pressure_torr);

    flux_Arp = flux_ions(armass, lambdai);
    flux_SiH3p = flux_ions(sih4mass, 2.9e-3);

    vsheath = 0.25*100.0;
}

void PlasmaChem::update_density(std::vector<double> pdens,
                                std::vector<double> ndens) {
    this->ndens = ndens;
    this->pdens = pdens;
}

int PlasmaChem::solve(const double ti,
			  const double tf,
			  const double dt,
			  const std::vector<double> &density_sourcedrain,
			  const double energy_sourcedrain,
			  const std::vector<double> &olddensity,
			  std::vector<double> &newdensity,
              double &min_dtq) {
    pdens = olddensity;
    ndens = olddensity;
    this->density_sourcedrain = density_sourcedrain;

    ssize_t dens_size = pdens.size();
    double atol[dens_size], rtol[dens_size];
    for(unsigned int q=0; q<dens_size; ++q) {
        atol[q] = 1.0e-6;
        rtol[q] = 1.0e-6;
    }

    // INTEGRATE LSODA begins
    // set options
    struct lsoda_opt_t plasma_opt = {0};
    //opt.h0 = ldtq;                // initial time step
    plasma_opt.ixpr = 0;                 // additional printing
    plasma_opt.rtol = rtol;              // relative tolerance
    plasma_opt.atol = atol;              // absolute tolerance
    plasma_opt.itask = 1;                // normal integration
    plasma_opt.mxstep = 100000000;      // max steps
    // set lsoda context
    struct lsoda_context_t plasma_ctx = {
                    .function = plasma_system,
                    .data = this,
                    .neq = static_cast<int>(dens_size),
                    .state = 1,
    };
    lsoda_prepare(&plasma_ctx, &plasma_opt);
    // time for lsoda
    double tin = ti;
    double tout = tf;
    // integrate
    lsoda(&plasma_ctx, &ndens[0], &tin, tout);
    min_dtq = std::min(min_dtq, plasma_ctx.common->hu);
    if (plasma_ctx.state <= 0) {
        std::cerr << "\nError in plasma solver\n";
        std::cerr << "\nerror istate = " << plasma_ctx.state;
        std::cerr << "\nmin dtq = " << min_dtq << "\n\n";
        std::terminate();
    }
    // free context
    lsoda_free(&plasma_ctx);

    newdensity = ndens;
    
    return 0;
}

double PlasmaChem::thermal_velocity(double mass) {
    return sqrt(2.0*KbTg/mass);
}

double PlasmaChem::diffusion_neutrals(double mass, double lambda) {
    return KbTg*lambda/(mass*thermal_velocity(mass));
}

double PlasmaChem::center2edge_neutrals(double mass) {
    double pf = 1.0 + (length/2)*thermal_velocity(mass)
              / (4.0*diffusion_neutrals(mass));
    return 1.0/pf;
}

double PlasmaChem::flux_neutrals(double mass) {
    return 0.25 * center2edge_neutrals(mass) * thermal_velocity(mass);
}

double PlasmaChem::bohm_velocity(double mass) {
    return sqrt(KbTg/mass);
}

double PlasmaChem::center2edge_ions(double mass, double lambda){
    double pf = sqrt(3.0+(0.5*length/lambda));
    return 1.0/pf;
}

double PlasmaChem::flux_ions(double mass, double lambda) {
    return center2edge_ions(mass, lambda) * bohm_velocity(mass);
}

double PlasmaChem::ion_velocity(double mass) {
    return sqrt(8.0*KbTg/(pi*mass));
}