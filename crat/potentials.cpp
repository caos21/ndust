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

  //////////////////////////////////////////////////////////////////////////////

  //--------------------------------------------------------------------------//

  //////////////////////////////////////////////////////////////////////////////



  
  // Compute potentials as function of ratios
  {
    // sizes
    unsigned int rsize = 100;
    unsigned int qsize = 302;

    // radius ratio
    darray r21array = linear<double>(0.5/50.0, 1.0, rsize);
    // charge ratio
    darray q21array = linear<double>(-5/296.0, 1.0, qsize);

    double eps = 11.68;

    for (unsigned int ir=0; ir<rsize; ++ir) {
      cout << '\t' << q21array[ir];
    }
    
    boost_array2d ipa_potentials;
    bgrid2d potgrid = {{rsize, qsize}};
    ipa_potentials.resize(potgrid);

    boost_array2d coulomb_potentials;
    coulomb_potentials.resize(potgrid);

    boost_array2d mpc_potentials;
    mpc_potentials.resize(potgrid);

    boost_array2d mpc_nterms;
    mpc_nterms.resize(potgrid);
    
    H5::H5File h5obj;
    int err = open_hdf5("/home/ben/git/ndust/data/potentials-vs-ratios.h5",
			h5obj, "create");
    if(err!=0) {
      BOOST_LOG_SEV(lg, error) << "h5 I/O error opening file. Terminate.";
      std::terminate();
    }

    string gname("Potentials");

    err = create_attrib_hdf5(h5obj,
			     gname,
			     "eps",
			     eps);
    
    err = create_group_hdf5(h5obj, gname);

    short nmin = 25;
    short nmax = 2000;
    short nstep = 5;

#pragma omp parallel for collapse(2)
    for (unsigned int ir=0; ir<rsize; ++ir) {
      for (unsigned int jq=0; jq<qsize; ++jq) {
	double r21 = r21array[ir];
	double q21 = q21array[jq];
	double rt = 1.0 + r21;
	//double rtmax = 2.5;//r1o - rt;
	//double q21 = q2o/q1o;
	
	potential_coulomb_funct coulombfunct(r21, q21);
	potential_ipa_funct pipafunct(r21, q21, eps);
	//potential_mpc_funct pmpcfunct(r21, q21, eps, 40);

	//double potential_prefactor = 1.0;//q1o*q1o/r1o;
	double pipa_rt = pipafunct(rt);
	ipa_potentials[ir][jq] = pipa_rt;

	double pcoulomb_rt = coulombfunct(rt);
	coulomb_potentials[ir][jq] = pcoulomb_rt;

	double final_pmpc_rt = 0.0;	
	double pcomp = pipa_rt;
	
	unsigned int initer=0;

	short nterms = 0;
	for(short n=nmin; n <= nmax; n+=nstep, ++initer) {
      
	  // mpc functor
	  potential_mpc_funct pmpcfunct(r21, q21, eps, n);

	  // mpc at contact
	  double pmpc_rt = pmpcfunct(rt);

	  double error_comp = max_pct_error(pcomp, pmpc_rt);

	  if ((error_comp < MPC_ERROR) && (initer>0)) {
	    final_pmpc_rt = pmpc_rt;
	    nterms = n;
	    n=nmax;
	  }
	  else {
	    if (n>nmax-nstep) {
	      std::cerr << "\n[ww] Max iterations exceeded\n";
	      final_pmpc_rt = pmpc_rt;
	      nterms = n;
	    }
	  }
	  pcomp = pmpc_rt;
	}

	mpc_potentials[ir][jq] = final_pmpc_rt;
	mpc_nterms[ir][jq] = nterms;
	//cout << '\t' << potentials[ir][jq];
      }
    }

    string sname_ipa("ipa_contact");
    err = write_dataset2d_hdf5(h5obj,
			       gname,
			       sname_ipa,
			       ipa_potentials,
			       rsize,
			       qsize);

    string sname_coulomb("coulomb_contact");
    err = write_dataset2d_hdf5(h5obj,
			       gname,
			       sname_coulomb,
			       coulomb_potentials,
			       rsize,
			       qsize);

    string sname_mpc("mpc_contact");
    err = write_dataset2d_hdf5(h5obj,
			       gname,
			       sname_mpc,
			       mpc_potentials,
			       rsize,
			       qsize);

    string sname_mpc_nterms("mpc_nterms");
    err = write_dataset2d_hdf5(h5obj,
			       gname,
			       sname_mpc_nterms,
			       mpc_nterms,
			       rsize,
			       qsize);

    string sname_r21("r21");
    err = write_dset_hdf5< darray >(h5obj,
				    gname,
				    sname_r21,
				    r21array,
				    rsize);

    string sname_q21("q21");
    err = write_dset_hdf5< darray >(h5obj,
				    gname,
				    sname_q21,
				    q21array,
				    qsize);
    

    //----------------------------------------------------------------------------
    //  IPA Barrier
    //----------------------------------------------------------------------------

    boost_array2d ipa_barrier;
    ipa_barrier.resize(potgrid);

    boost_array2d ipa_rbarrier;
    ipa_rbarrier.resize(potgrid);

    boost_array2d mpc_barrier_hybrid;
    mpc_barrier_hybrid.resize(potgrid);    

#pragma omp parallel for collapse(2)
    for (unsigned int ir=0; ir<rsize; ++ir) {
      for (unsigned int jq=0; jq<qsize; ++jq) {
	double r21 = r21array[ir];
	double q21 = q21array[jq];
	double rt = 1.0 + r21;

	double rmin = rt;
	double rmax = 100.0*rt;// 100 times contact radii

	double min = rt;
	double max = rmax;

	// nterms for hybrid approach
	short nterms = mpc_nterms[ir][jq];

	force_ipa_funct fipafunct(r21, q21, eps);
	// force ipa at contact 
	double fipa_rmin = fipafunct(rmin);
	// Force at r max
	double fipa_rmax = fipafunct(rmax);

	// checks if minimum exists
	//	#pragma omp ordered
	if (fipa_rmin*fipa_rmax < 0.0) {	  
	  boost::uintmax_t bmax_iter = ROOT_MAXITER;
	  tools::eps_tolerance<double> tol = ROOT_TOL;

	  std::pair<double, double> pair_fipa;
	  try {
	    pair_fipa = tools::toms748_solve(fipafunct, min, max, fipa_rmin, fipa_rmax, tol, bmax_iter);
	  }
	  catch(const std::exception& exc) {
	    std::cerr << '\n' << exc.what() << '\n';
	    //std::terminate();
	    std::cerr << "\n[ee] No barrier\n";
	  }
	  if(bmax_iter > 990){
	    std::cerr << "\n ERROR max iter " << bmax_iter << "\n\n";
	    std::terminate();
	  }
	  double rbarrier = 0.5*(pair_fipa.first+pair_fipa.second);
	  if(rbarrier>=min){
	    //*********************** USE SAME COEFFICIENTS OF FORCE
	    // ipa functor
	    potential_ipa_funct pipafunct(r21, q21, eps);      
	    // ipa at barrier
	    double pipa_barrier = pipafunct(rbarrier);
	    // mpc functor
	    potential_mpc_funct pmpcfunct(r21, q21, eps, nterms);
	    // mpc at ipa rbarrier
	    double pmpc_barrier = pmpcfunct(rbarrier);
#pragma omp critical
	    {
	      ipa_barrier[ir][jq] = pipa_barrier;
	      ipa_rbarrier[ir][jq] = rbarrier;
	      
	      // mpc at contact
	      mpc_barrier_hybrid[ir][jq] = pmpc_barrier;
	    }
	  }
	  else{
	    std::cerr << "\n ERROR Negative rbarrier " << rbarrier << '\n';
	    std::terminate();
	  }	
	} 
      }
    }
    

    string sname_ipa_barrier("ipa_barrier");
    err = write_dataset2d_hdf5(h5obj,
			       gname,
			       sname_ipa_barrier,
			       ipa_barrier,
			       rsize,
			       qsize);

    string sname_ipa_rbarrier("ipa_rbarrier");
    err = write_dataset2d_hdf5(h5obj,
			       gname,
			       sname_ipa_rbarrier,
			       ipa_rbarrier,
			       rsize,
			       qsize);

    string sname_mpc_barrier_hybrid("mpc_barrier_hybrid");
    err = write_dataset2d_hdf5(h5obj,
			       gname,
			       sname_mpc_barrier_hybrid,
			       mpc_barrier_hybrid,
			       rsize,
			       qsize);

    //----------------------------------------------------------------------------
    //  MPC Barrier
    //----------------------------------------------------------------------------

    boost_array2d mpc_barrier;
    mpc_barrier.resize(potgrid);

    boost_array2d mpc_rbarrier;
    mpc_rbarrier.resize(potgrid);

#pragma omp parallel for collapse(2) schedule(nonmonotonic:dynamic)
    for (unsigned int ir=0; ir<rsize; ++ir) {
      for (unsigned int jq=0; jq<qsize; ++jq) {
	double r21 = r21array[ir];
	double q21 = q21array[jq];
	double rt = 1.0 + r21;

	double rmin = rt;
	double rmax = 100.0*rt;// 100 times contact radii

	double min = rt;
	double max = rmax;

	short nterms = mpc_nterms[ir][jq];
	
	force_mpc_funct fmpcfunct(r21, q21, eps, nterms);

	// force mpc at contact 
	double fmpc_rmin = fmpcfunct(rmin);
	// Force at r max
	double fmpc_rmax = fmpcfunct(rmax);

	// checks if minimum exists
	if (fmpc_rmin*fmpc_rmax < 0.0) {
	  // std::cerr << "\n[ii] Mixed phi_rt = " << fmpc_rmin << '\t' << fmpc_rmax;
	  boost::uintmax_t bmax_iter = ROOT_MAXITER;
	  tools::eps_tolerance<double> tol = ROOT_TOL;

	  std::pair<double, double> pair_fmpc;
	  // try {
	  pair_fmpc = tools::toms748_solve(fmpcfunct, min, max, fmpc_rmin, fmpc_rmax, tol, bmax_iter);
	  if(bmax_iter > 990){
	    std::cerr << "\n ERROR max iter " << bmax_iter << "\n\n";
	    std::terminate();
	  }
	  double rbarrier = 0.5*(pair_fmpc.first+pair_fmpc.second);
	  if(rbarrier>=min){
	    //*********************** USE SAME COEFFICIENTS OF FORCE
#pragma omp critical
	    {
	      // mpc functor
	      potential_mpc_funct pmpcfunct(r21, q21, eps, nterms);
	      // mpc at contact
	      double pmpc_barrier = pmpcfunct(rbarrier);
	      mpc_barrier[ir][jq] = pmpc_barrier;
	      mpc_rbarrier[ir][jq] = rbarrier;
	    }
	  }
	  else{
	    std::cerr << "\n ERROR Negative rbarrier " << rbarrier << '\n';
	    std::terminate();
	  }	
	}	
      }
    }

    string sname_mpc_barrier("mpc_barrier");
    err = write_dataset2d_hdf5(h5obj,
			       gname,
			       sname_mpc_barrier,
			       mpc_barrier,
			       rsize,
			       qsize);

    string sname_mpc_rbarrier("mpc_rbarrier");
    err = write_dataset2d_hdf5(h5obj,
			       gname,
			       sname_mpc_rbarrier,
			       mpc_rbarrier,
			       rsize,
			       qsize);


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

