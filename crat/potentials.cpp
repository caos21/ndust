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



int main(int argc, char* argv[], char* envp[]) {

  Logger tlog("potentials");
  
  src::severity_logger< severity_level > lg = tlog.get_lg();  

  {
    // particle 1 radii and charge
    unsigned int rtsize = 100;
    unsigned int r1size = 9;
    darray r1array = linear<double>(3.5e-9, 4.5e-9, r1size);
    double q1o = -eCharge;
    double r2o = 2.25e-9;
    double q2o = -3*eCharge;
    double eps = 11.68;

    
    boost_array2d potentials;    
    bgrid2d potgrid = {{rtsize, 4}};
    potentials.resize(potgrid);

    H5::H5File h5obj;
    int err = open_hdf5("potentials.h5", h5obj, "create");
    if(err!=0) {
      BOOST_LOG_SEV(lg, error) << "h5 I/O error opening file. Terminate.";
      std::terminate();
    }

    string gname("Potentials");

    err = create_attrib_hdf5(h5obj,
			     gname,
			     "q1",
			     q1o);
    
    err = create_attrib_hdf5(h5obj,
			     gname,
			     "r2",
			     r2o);

    err = create_attrib_hdf5(h5obj,
			     gname,
			     "q2",
			     q2o);

    err = create_attrib_hdf5(h5obj,
			     gname,
			     "eps",
			     eps);
    
    err = create_group_hdf5(h5obj, gname);
    
    for (unsigned int i=0; i<r1size; ++i) {
      double r1o = r1array[i];
      double r21 = r2o/r1o;
      double rt = 1 + r21;
      double rtmax = 2.5;//r1o - rt;
      double q21 = q2o/q1o;
      cout << '\n' << r1array[i]
	   << '\t' << r21
	   << '\t' << q21
	   << '\t' << rt
	   << '\t' << rtmax;
	
      darray rtarray = linear<double>(rt, rtmax, rtsize);

      potential_coulomb_funct coulombfunct(r21, q21);
      potential_ipa_funct pipafunct(r21, q21, eps);
      potential_mpc_funct pmpcfunct(r21, q21, eps, 40);
      double potential_prefactor = q1o*q1o/r1o;
      
      for (unsigned int j=0; j<rtsize; ++j) {
	double rd = rtarray[j];
	potentials[j][0] = rd*r1o;
	potentials[j][1] = potential_prefactor * coulombfunct(rd);
	potentials[j][2] = potential_prefactor * pipafunct(rd);
	potentials[j][3] = potential_prefactor * pmpcfunct(rd);

      }

      string sname = to_string(i);
      err = write_dataset2d_hdf5(h5obj,
                                 gname,
                                 sname,
                                 potentials,
                                 rtsize,
                                 4);

      err = create_dsattrib_hdf5<double>(h5obj,
					 gname,
					 sname,
					 "r1",
					 r1o);
      
      cout << '\n' << gname << '\t' << sname;      


    }

    hid_t id = h5obj.getId();
    
    err = close_hdf5(h5obj);
    if(err!=0) {
      BOOST_LOG_SEV(lg, error) << "h5 I/O error closing file. Terminate.";
      std::terminate();
    }
    BOOST_LOG_SEV(lg, info) << "Closed file Id: " << id;

  }
  
  return 0;
}

