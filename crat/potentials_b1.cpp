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

#define PTOL 1e-6

#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <functional>
#include <cstdlib>

#include <omp.h>

#include "io.h"
#include "../include/log.h"
#include "../include/array.h"
#include "../include/h5plasma.h"
#include "../include/eint.h"
#include "../include/constants.h"
#include "../include/enhancement.h"

using namespace std;
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



int main(int argc, char* argv[], char* envp[]) {

  Logger tlog("potentials");
  
  src::severity_logger< severity_level > lg = tlog.get_lg();  

  // get num threads
  int num_threads = omp_get_num_threads();
  BOOST_LOG_SEV(lg, info) << "Number of threads = " << num_threads;
  
  // get max threads
  int max_num_threads = omp_get_max_threads();
  BOOST_LOG_SEV(lg, info) <<  "Max number of threads = " << max_num_threads;

  char *input_opt;
  input_opt = get_cmd_option(argv, argv + argc, "-t");
  if(input_opt) {
    // get num threads from input
    num_threads = std::stoi(input_opt);
    // checks if num_threads > max_num_threads
    if(num_threads>max_num_threads) {
      num_threads = 1;
    }
    // set num_threads
    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
      #pragma omp single
      BOOST_LOG_SEV(lg, info) << "Set number of threads = " << num_threads;
    }
  }

  BOOST_LOG_SEV(lg, info) <<  "Get omp nested = " << omp_get_nested();
  BOOST_LOG_SEV(lg, info) <<  "Get omp cancellation = " << omp_get_cancellation();

  //
  
  // Compute potentials for differents r1
  {
    // particle 1 radii and charge
    double r1o = 50.0e-9;
    double q1o = -eCharge;
    double r2o = 1e-9;
    double q2o = -eCharge;
    double eps = 11.68;

    //double rto = r1o + r2o;
    
    unsigned terms = 10;
    unsigned tsize = 100;
    double r21 = r2o/r1o;
    double q21 = q2o/q1o;
    double rt = 1 + r21;    

    double potential_prefactor = (q1o*q1o/r1o)/eCharge;

    potential_ipa_funct pipafunct(r21, q21, eps);
    potential_mpc_funct pmpcfunct(r21, q21, eps, terms);
      
    double prev_val = potential_prefactor * pmpcfunct(rt);
    double new_val = potential_prefactor * pmpcfunct(rt);
    double error = 0.0;
    
    for (unsigned int i=0; i<tsize; ++i) {

      terms += 10;	
      potential_mpc_funct pmpcfunct(r21, q21, eps, terms);

      new_val = potential_prefactor * pmpcfunct(rt);

      error = 100.0* fabs((new_val-prev_val)/new_val);

      cerr << '\n' << terms
	   << '\t' << potential_prefactor * pipafunct(rt)
	   << '\t' << new_val
	   << '\t' << error;
	  
      prev_val = new_val;

      if(error<0.1) break;

    }

    unsigned rtsize = 100;
    darray rtarray = linear<double>(rt, rt*4.0, 100);

    potential_mpc_funct pmpcfunct2(r21, q21, eps, terms);
	  
    for (unsigned int j=0; j<rtsize; ++j) {
        double rd = rtarray[j];
	double rdnm = rd*r1o;
	//potentials[j][1] = potential_prefactor * coulombfunct(rd);
	//potentials[j][2] = potential_prefactor * pipafunct(rd);
	//potentials[j][3] = potential_prefactor * pmpcfunct(rd);
	cout << '\n' << rdnm
	     << '\t' << potential_prefactor * pipafunct(rd)
	     << '\t' << potential_prefactor * pmpcfunct2(rd);
      }
  }

  cout << '\n';
}

