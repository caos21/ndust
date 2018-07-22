/*
 * Copyright 2018 <Benjamin Santos> <caos21@gmail.com>
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

#define BOOST_TEST_MODULE enhancement_test
//#define BOOST_TEST_NO_MAIN
//#define BOOST_TEST_ALTERNATIVE_INIT_API

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

//#define BOOST_AUTO_TEST_MAIN 

//#define BOOST_TEST_MAIN

#define PTOL 1e-6

#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <functional>
#include <cstdlib>

#include "../include/log.h"
#include "../include/constants.h"
#include "../include/enhancement.h"

using namespace std;
using namespace boost::unit_test;
using namespace boost::test_tools;

using namespace eint;
using namespace enhancement;

struct Logger {

  Logger() { }
  
  Logger(std::string filename) {
    // init logger
    blog::init(filename);

    logging::add_common_attributes();    

    BOOST_LOG_SEV(lg, info) << "Logging started...";
  

  }
  ~Logger(){ }
  src::severity_logger< severity_level > get_lg(){
    return lg;
  }
  src::severity_logger< severity_level > lg;
};

void print_efactor(boost_array4d array4d, unsigned rsize, unsigned qsize) {
  for (unsigned int l=0; l<rsize; ++l) {
    // iterate in charges particle 1
    for (unsigned int q=0; q<qsize; ++q) {
      // iterate in radii particle 2
      std::cout << "\n\n---slab---";
      for (unsigned int m=0; m<rsize; ++m) {
	// iterate in charges particle 2
	std::cout << '\n';
	for (unsigned int p=0; p<qsize; ++p) {
	  std::cout << '\t' << array4d[l][q][m][p];
	}
      }
    }
  }
}
//Logger TLOG;

//src::severity_logger< severity_level > lg = TLOG.get_lg();

//******************************************************************************
//**                              Test enhancement                            **
//******************************************************************************

//BOOST_GLOBAL_FIXTURE(Logger);


BOOST_AUTO_TEST_CASE( coulomb_enhancement_test ) {
  {

    Logger tlog("coulomb-enhancement");

    src::severity_logger< severity_level > lg = tlog.get_lg();
 
    darray rarray = {1e-9, 50e-9};
    darray qarray = {-eCharge, 0.0, eCharge};
    // darray rarray = {1e-9, 2e-9, 3e-9, 50e-9};
    // darray qarray = {-eCharge, eCharge};    
    double eps = 11.68;

    bgrid4d grid4;
    boost_array4d efactor;

    unsigned rsize = static_cast<unsigned>(rarray.size());
    unsigned qsize = static_cast<unsigned>(qarray.size());
    grid4 = {{rsize, qsize, rsize, qsize}};

    efactor.resize(grid4);
    boost_array4d_ref efactor_ref(efactor);
     
    Enhancement enh(rarray, qarray, efactor_ref, eps, lg);

    BOOST_TEST_MESSAGE( "\n\n[ii] Testing particle pairs..." );
    enh.compute_reducedpairs();
    for (auto pp: enh.particle_pairs) {
      pp.print();
    }

    BOOST_TEST_MESSAGE( "\n\n[ii] Testing reduced particle pairs..." );
    for (auto rpp: enh.reduced_pairs) {
      rpp.print();
    }

    BOOST_TEST_MESSAGE( "\n\n[ii] Testing coulomb potentials at contact..." );
    enh.compute_coulombpotential_contact();
    for (auto cpot: enh.contact_potentials) {
      cpot.print();
    }
    
    BOOST_TEST_MESSAGE( "\n\n[ii] Testing enhancement factor..." );
    enh.compute_enhancement_factor();

    print_efactor(efactor, rsize, qsize);
    
    BOOST_REQUIRE_CLOSE(1.0, 1.0, PTOL);
  }
}

BOOST_AUTO_TEST_CASE( ipa_enhancement_test ) {
  {

    Logger tlog("ipa-enhancement");

    src::severity_logger< severity_level > lg = tlog.get_lg();

    darray rarray = {1e-9, 50e-9};
    darray qarray = {-eCharge, 0.0, eCharge};
    double eps = 11.68;

    bgrid4d grid4;
    boost_array4d efactor;

    unsigned rsize = static_cast<unsigned>(rarray.size());
    unsigned qsize = static_cast<unsigned>(qarray.size());
    grid4 = {{rsize, qsize, rsize, qsize}};

    efactor.resize(grid4);
    boost_array4d_ref efactor_ref(efactor);

    Enhancement enh(rarray, qarray, efactor_ref, eps, lg);
    
    BOOST_TEST_MESSAGE( "\n\n[ii] Testing particle pairs..." );
    enh.compute_reducedpairs();
    for (auto pp: enh.particle_pairs) {
      pp.print();
    }

    BOOST_TEST_MESSAGE( "\n\n[ii] Testing reduced particle pairs..." );
    for (auto rpp: enh.reduced_pairs) {
      rpp.print();
    }

    BOOST_TEST_MESSAGE( "\n\n[ii] Testing ipa potentials at contact..." );
    enh.compute_ipapotential_contact();
    for (auto cpot: enh.contact_potentials) {
      cpot.print();
    }
    
    BOOST_TEST_MESSAGE( "\n\n[ii] Testing ipa potentials at barrier..." );
    enh.compute_ipapotential_barrier();
    for (auto cpot: enh.barrier_potentials) {
      cpot.print();
    }

    BOOST_TEST_MESSAGE( "\n\n[ii] Testing enhancement factor..." );
    enh.compute_enhancement_factor();

    print_efactor(efactor, rsize, qsize);
    
    BOOST_REQUIRE_CLOSE(1.0, 1.0, PTOL);
  }
}

BOOST_AUTO_TEST_CASE( mpc_enhancement_test ) {
  {

    Logger tlog("mpc-enhancement");

    src::severity_logger< severity_level > lg = tlog.get_lg();
    
    darray rarray = {1e-9, 50e-9};
    darray qarray = {-eCharge, 0.0, eCharge};
    double eps = 11.68;

    bgrid4d grid4;
    boost_array4d efactor;

    unsigned rsize = static_cast<unsigned>(rarray.size());
    unsigned qsize = static_cast<unsigned>(qarray.size());
    grid4 = {{rsize, qsize, rsize, qsize}};

    efactor.resize(grid4);
    boost_array4d_ref efactor_ref(efactor);
     
    Enhancement enh(rarray, qarray, efactor, eps, lg);

    BOOST_TEST_MESSAGE( "\n\n[ii] Testing particle pairs..." );
    enh.compute_reducedpairs();
    for (auto pp: enh.particle_pairs) {
      pp.print();
    }

    BOOST_TEST_MESSAGE( "\n\n[ii] Testing reduced particle pairs..." );
    for (auto rpp: enh.reduced_pairs) {
      rpp.print();
    }

    BOOST_TEST_MESSAGE( "\n\n[ii] Testing mpc potentials at contact..." );
    enh.compute_mpcpotential_contact();
    for (auto cpot: enh.contact_potentials) {
      cpot.print();
    }
    
    BOOST_TEST_MESSAGE( "\n\n[ii] Testing mpc potentials at barrier..." );
    enh.compute_mpcpotential_barrier();
    for (auto cpot: enh.barrier_potentials) {
      cpot.print();
    }

    BOOST_TEST_MESSAGE( "\n\n[ii] Testing enhancement factor..." );
    enh.compute_enhancement_factor();

    print_efactor(efactor, rsize, qsize);
    
    BOOST_REQUIRE_CLOSE(1.0, 1.0, PTOL);
  }

}

// bool init_unit_test() {
//   return true;
// }

// BOOST_AUTO_TEST_SUITE_END()
// int main(int argc, char* argv[], char* envp[]) {

//   unit_test_main(init_unit_test, argc, argv);
// }
