/*
 * Copyright 2017 <Benjamin Santos> <caos21@gmail.com>
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

#include "NEvo.h"

int NEvo::open() {
  std::string fullname = dirname + grid_filename + ".h5";

  int err = open_hdf5(fullname, h5obj, "write");
  if(err!=0) {
    BOOST_LOG_SEV(lg, error) << "h5 I/O error opening file. Terminate.";
    std::terminate();
  }
  BOOST_LOG_SEV(lg, info) << "Open file for read/write: "
                          << fullname << " -> Success. Id: "
                          << h5obj.getId();

  fullname = dirname + plasma_filename + ".h5";
  err = open_hdf5(fullname, h5obj_plasma, "write");
  if(err!=0) {
    BOOST_LOG_SEV(lg, error) << "h5 I/O error opening file. Terminate.";
    std::terminate();
  }
  BOOST_LOG_SEV(lg, info) << "Open file for read/write: "
                          << fullname << " -> Success. Id: "
                          << h5obj_plasma.getId();

  fullname = dirname + nano_filename + ".h5";
  err = open_hdf5(fullname, h5obj_nano, "write");
  if(err!=0) {
    BOOST_LOG_SEV(lg, error) << "h5 I/O error opening file. Terminate.";
    std::terminate();
  }
  BOOST_LOG_SEV(lg, info) << "Open file for read/write: "
                          << fullname << " -> Success. Id: "
                          << h5obj_nano.getId();

  return 0;
}

int NEvo::close() {
  delete(qsys);
  
  delete(ls_qsys);
  
  //delete(sol);

  hid_t id = h5obj.getId();

  int err = close_hdf5(h5obj);
  if(err!=0) {
    BOOST_LOG_SEV(lg, error) << "h5 I/O error closing file. Terminate.";
    std::terminate();
  }
  BOOST_LOG_SEV(lg, info) << "Closed file Id: " << id;

  id = h5obj_plasma.getId();
  err = close_hdf5(h5obj_plasma);
  if(err!=0) {
    BOOST_LOG_SEV(lg, error) << "h5 I/O error closing file. Terminate.";
    std::terminate();
  }
  BOOST_LOG_SEV(lg, info) << "Closed file Id: " << id;

  id = h5obj_nano.getId();
  err = close_hdf5(h5obj_nano);
  if(err!=0) {
    BOOST_LOG_SEV(lg, error) << "h5 I/O error closing file. Terminate.";
    std::terminate();
  }
  BOOST_LOG_SEV(lg, info) << "Closed file Id: " << id;

  moments_file->close();
  delete(moments_file);
  BOOST_LOG_SEV(lg, info) << "Closed file moments";

  if(sih4_file != NULL) {
    sih4_file->close();
  }
  delete(sih4_file);
  BOOST_LOG_SEV(lg, info) << "Closed file sih4";
  
  return 0;
}

int NEvo::read() {
  hid_t id = h5obj.getId();

  cr = CRate(h5obj, lg);
  cr.read_results();

  pm = PlasmaModel(h5obj_plasma, lg);
  pm.read();

  nm = NanoModel(h5obj_nano, lg);
  nm.read();

  return 0;
}

int NEvo::write() {
  hid_t id = h5obj.getId();

//   BOOST_LOG_SEV(lg, info) << "Writing eta creation factor";
//   write_etafactor();
// 
//   BOOST_LOG_SEV(lg, info) << "Writing death factor";
//   write_deathfactor();
// 
//   BOOST_LOG_SEV(lg, info) << "Writing enhancement factor";
//   write_efactor();
// 
//   BOOST_LOG_SEV(lg, info) << "Writing coagulation rate";
//   write_rcoagulation();

  return 0;
}

int NEvo::evolve() {
  // IMPORTANT set nested parallelism
  omp_set_nested(true);
#pragma omp single
  {// begin master omp, controls the iterations in t (time)
    
    BOOST_LOG_SEV(lg, info) << "Get number of threads = " << omp_get_num_threads();
    BOOST_LOG_SEV(lg, info) << "Nested = " << omp_get_nested();
    BOOST_LOG_SEV(lg, info) << "Dynamic = " << omp_get_dynamic();
  }
  
//   check = 2;
//   // initial conditions, in boost vector_type x( 2 , 1.0 );
//   N_Vector xini = NULL;
//   xini = N_VNew_Serial(cr.gm.chrgs.nsections);
//   NV_Ith_S(xini, 0) = RCONST(2.0);
//   NV_Ith_S(xini, 1) = RCONST(1.0);

//     std::cerr << "\n Initial condition : "
//               << NV_Ith_S(xini,0) << '\t' << NV_Ith_S(xini,1) << '\n' << NV_LENGTH_S(xini) << "\n"; 

//   // Instantiate the solver
//   sol = new Solver(this, xini, 0.0);
// //  Solver sol(xini, 0.0);
//   check = 3;
//   sol->print();
//   // final time
//   realtype tf = RCONST(50.0);
//   // delta time
//   realtype dt = RCONST(0.01);
// //   sol->compute(tf, dt);
// //   sol.compute(tf, dt);


//   // or solve one step
// //   sol.compute_step(0.0, dt);

//   // free memory for xini
//   N_VDestroy_Serial(xini);

  int err;

  unsigned int npcount = 0;

  err = compute_precompute();

  qsys = new qsystem(this);
  
  auto stepper = odeint::make_controlled( 1.0e-6, 1.0e-6,
                                          odeint::runge_kutta_dopri5< state_type >() );


  clock_t begin_sim = std::clock();
  ctime = 0.0;
  for(unsigned int t=0; ; ++t, ++npcount) {
    ctime += nm.tm.ndeltat;

    clock_t one_step_begin = std::clock();

    evolve_one_step(ctime);

    if(npcount>static_cast<unsigned int>((0.01*nm.tm.tstop)/nm.tm.ndeltat)) {
      clock_t elapsed_sim = std::clock();
      double elapsed_secs = double(elapsed_sim - begin_sim) / CLOCKS_PER_SEC;
      BOOST_LOG_SEV(lg, info) << "Writing results at simulation time = "
                              << ctime << " [elapsed secs = " 
                              << elapsed_secs << "]";
      err=write_partial_results(ctime);
      clock_t one_step_end = std::clock();
      double elapsed_secs_onestep  = double(one_step_end - one_step_begin) / CLOCKS_PER_SEC;
      BOOST_LOG_SEV(lg, info) << "One step elapsed secs = "  << elapsed_secs_onestep;
      npcount=0;
    }

    if(ctime>nm.tm.tstop-nm.tm.ndeltat) {
      BOOST_LOG_SEV(lg, info) << "Done!, time = " << ctime;
      clock_t end_sim = std::clock();
      double elapsed_secs = double(end_sim - begin_sim) / CLOCKS_PER_SEC;
      BOOST_LOG_SEV(lg, info) << "Nanoparticle module elapsed time "
                              << elapsed_secs;
      err=write_partial_results(ctime);
      err = write_dataset2d_hdf5(h5obj_nano, "Density", "density", ndens,
                                 cr.gm.vols.nsections, cr.gm.chrgs.nsections);

      return 0;
    }
    pdens = ndens;
  } 

  return -1;
}

// -------------------- evolve adaptive -----------------

int NEvo::evolve_radapt() {

  int err;

  // store the original nano delta time
  double orig_dtn = nm.tm.ndeltat;
  
  if(nm.rs.wsih4==1) {
    BOOST_LOG_SEV(lg, info) << "With SiH4 rates";
    err = compute_precompute_sih4();
  }
  else {
    err = compute_precompute();
  }
  
  // lsoda
  ls_qsys = new ls_qsystem(this);
  
  // IMPORTANT set nested parallelism
  omp_set_nested(true);
  
  clock_t begin_sim = std::clock();
  ctime = 0.0;
  err = -1;
  int err_onestep = 0;
  unsigned int write_count = static_cast<unsigned int>((0.01*nm.tm.tstop)
						       /orig_dtn);  
#pragma omp parallel //shared(ctime, nm, err)
  {//begin parallel omp
#pragma omp master
    {// begin master omp, controls the iterations in t (time)

      omp_info();
      unsigned int npcount = 0;
      // get num_threads
      int num_threads = omp_get_num_threads();
      // execute loop in master thread
      
      for(unsigned int t=0; err<0; ++t, ++npcount) {
	
	//clock_t one_step_begin = std::clock();
	double one_step_begin = omp_get_wtime();

#pragma omp parallel /*num_threads(num_threads) shared(ctime)*/
	{
	  err_onestep = evolve_one_step_adapt(ctime);
	}
	if(npcount>write_count) {
	  write_one_step(ctime, begin_sim, one_step_begin);
	  npcount=0;
	}
	if(ctime>nm.tm.tstop-nm.tm.ndeltat) {
	  BOOST_LOG_SEV(lg, info) << "Done!, time = " << ctime;
	  clock_t end_sim = std::clock();
	  double elapsed_secs = double(end_sim - begin_sim) / CLOCKS_PER_SEC;
	  BOOST_LOG_SEV(lg, info) << "Nanoparticle module elapsed time "
				  << elapsed_secs;
	  int err_ = write_partial_results(ctime);
	  err_ = write_dataset2d_hdf5(h5obj_nano, "Density", "density", ndens,
				      cr.gm.vols.nsections, cr.gm.chrgs.nsections);
#pragma omp critical
	  {
	    err = 0;
	  }	
	  //break;
	}
	if (err_onestep == 0) {
	  // update time
	  ctime += nm.tm.ndeltat;
	  // update density
	  pdens = ndens;
	  // WARNING reset dtn
	  //nm.tm.ndeltat = orig_dtn;
	  //BOOST_LOG_SEV(lg, info) << "Reset to original dtn = " << nm.tm.ndeltat;	    
	}
      }
    }// end master omp
  }// end parallel omp    

  return err;
}

//
// -------------------- evolve adaptive -----------------
//

// -------------------- evolve openmp --------------------

int NEvo::evolve_omp() {

  int err;

  // store the original nano delta time
  double orig_dtn = nm.tm.ndeltat;

  // // time delta is same for all sections
  // for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
  //   for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
  //     dt2d[l][q] = nm.tm.ndeltat;
  //   }
  // }
  
  if(nm.rs.wsih4==1) {
    BOOST_LOG_SEV(lg, info) << "With SiH4 rates";
    err = compute_precompute_sih4();

    // lsoda
    ls_qsys = new ls_qsystem(this);

    // IMPORTANT set nested parallelism
    omp_set_nested(true);
  
    clock_t begin_sim = std::clock();
    ctime = 0.0;
    err = -1;
    int err_onestep = 0;
#pragma omp parallel //shared(ctime, nm, err)
    {//begin parallel omp
#pragma omp master
      {// begin master omp, controls the iterations in t (time)
	BOOST_LOG_SEV(lg, info) << "Master thread";
	omp_info();
	
	unsigned int npcount = 0;
	// get num_threads
	int num_threads = omp_get_num_threads();
	// execute loop in master thread

	for(unsigned int t=0; err<0; ++t, ++npcount) {

	  //clock_t one_step_begin = std::clock();
	  double one_step_begin = omp_get_wtime();
// fork
#pragma omp parallel /*num_threads(num_threads) shared(ctime)*/
	  {
	    err_onestep = evolve_one_step_omp(ctime);
	  }
	  if(npcount>static_cast<unsigned int>((0.01*nm.tm.tstop)/orig_dtn)) {
	    write_one_step(ctime, begin_sim, one_step_begin);
	    npcount=0;
	  }
	  if(ctime>nm.tm.tstop-nm.tm.ndeltat) {
	    BOOST_LOG_SEV(lg, info) << "Done!, time = " << ctime;
	    clock_t end_sim = std::clock();
	    double elapsed_secs = double(end_sim - begin_sim) / CLOCKS_PER_SEC;
	    BOOST_LOG_SEV(lg, info) << "Nanoparticle module elapsed time "
				    << elapsed_secs;
	    int err_ = write_partial_results(ctime);
	    err_ = write_dataset2d_hdf5(h5obj_nano, "Density", "density", ndens,
					cr.gm.vols.nsections, cr.gm.chrgs.nsections);
#pragma omp critical
	    {
	      err = 0;
	    }	
	    //break;
	  }
	  if (err_onestep == 0) {
	    // update time
	    ctime += nm.tm.ndeltat;
	    // update density
	    pdens = ndens;
	    // WARNING reset dtn
	    //nm.tm.ndeltat = orig_dtn;
	    //BOOST_LOG_SEV(lg, info) << "Reset to original dtn = " << nm.tm.ndeltat;	    
	  }
	}
      }// end master omp
    }// end parallel omp    
  }
  else {
    // check = 2;
    // // initial conditions, in boost vector_type x( 2 , 1.0 );
    // N_Vector xini = NULL;
    // xini = N_VNew_Serial(cr.gm.chrgs.nsections);
    // NV_Ith_S(xini, 0) = RCONST(2.0);
    // NV_Ith_S(xini, 1) = RCONST(1.0);

    // std::cerr << "\n Initial condition : "
    //           << NV_Ith_S(xini,0) << '\t' << NV_Ith_S(xini,1) << '\n' << NV_LENGTH_S(xini) << "\n"; 

    // // Instantiate the solver
    // sol = new Solver(this, xini, 0.0);
    // //  Solver sol(xini, 0.0);
    // check = 3;
    // sol->print();
    // // final time
    // realtype tf = RCONST(50.0);
    // // delta time
    // realtype dt = RCONST(0.01);
    // //   sol->compute(tf, dt);
    // //   sol.compute(tf, dt);


    // // or solve one step
    // //   sol.compute_step(0.0, dt);

    // // free memory for xini
    // N_VDestroy_Serial(xini);

    err = compute_precompute();

    qsys = new qsystem(this);

    auto stepper = odeint::make_controlled( 1.0e-6, 1.0e-6,
					    odeint::runge_kutta_dopri5< state_type >() );

  
    // lsoda
    ls_qsys = new ls_qsystem(this);

    // IMPORTANT set nested parallelism
    omp_set_nested(true);
  
    clock_t begin_sim = std::clock();
    ctime = 0.0;
    err = -1;
#pragma omp parallel //shared(ctime, nm, err)
    {//begin parallel omp
#pragma omp master
      {// begin master omp, controls the iterations in t (time)

	BOOST_LOG_SEV(lg, info) << "Set number of threads = " << omp_get_num_threads();
	BOOST_LOG_SEV(lg, info) << "Nested = " << omp_get_nested();
	BOOST_LOG_SEV(lg, info) << "Dynamic = " << omp_get_dynamic();
	BOOST_LOG_SEV(lg, info) << "Max active levels " << omp_get_max_active_levels();
	unsigned int npcount = 0;
	// get num_threads
	int num_threads = omp_get_num_threads();
	// execute loop in master thread
	double writing_step_begin = 0.0;
	for(unsigned int t=0; err<0; ++t, ++npcount) {
	  // update time
	  ctime += nm.tm.ndeltat;

	  //clock_t one_step_begin = std::clock();
	  double one_step_begin = omp_get_wtime();

	  //BOOST_LOG_SEV(lg, info) << "Number threads " << omp_get_num_threads();
#pragma omp parallel /*num_threads(num_threads) shared(ctime)*/
	  {
	    evolve_one_step_omp(ctime);
	  }
	  if(npcount>static_cast<unsigned int>((0.01*nm.tm.tstop)/nm.tm.ndeltat)) {
	    clock_t elapsed_sim = std::clock();
	    double elapsed_secs = double(elapsed_sim - begin_sim) / CLOCKS_PER_SEC;
	    BOOST_LOG_SEV(lg, info) << "Writing results at simulation time = "
				    << ctime << " [elapsed secs = " 
				    << elapsed_secs << "]";
	    int err_ = write_partial_results(ctime);
	    //clock_t one_step_end = std::clock();
	    double one_step_end = omp_get_wtime();
	    //double elapsed_secs_onestep  = double(one_step_end - one_step_begin) / CLOCKS_PER_SEC;
	    double elapsed_secs_onestep = one_step_end - one_step_begin;
	    BOOST_LOG_SEV(lg, info) << "One step elapsed secs = "  << elapsed_secs_onestep;

	    double writing_step_end = omp_get_wtime();
	    double elapsed_secs_writing = writing_step_end - writing_step_begin;
	    BOOST_LOG_SEV(lg, info) << "Writing time/Delta elapsed sec = " << ctime 
				    << "\t" << elapsed_secs_writing;    
	    // reset clock
	    writing_step_begin = omp_get_wtime();
	    npcount=0;
	  }
	  if(ctime>nm.tm.tstop-nm.tm.ndeltat) {
	    BOOST_LOG_SEV(lg, info) << "Done!, time = " << ctime;
	    clock_t end_sim = std::clock();
	    double elapsed_secs = double(end_sim - begin_sim) / CLOCKS_PER_SEC;
	    BOOST_LOG_SEV(lg, info) << "Nanoparticle module elapsed time "
				    << elapsed_secs;
	    int err_ = write_partial_results(ctime);
	    err_ = write_dataset2d_hdf5(h5obj_nano, "Density", "density", ndens,
					cr.gm.vols.nsections, cr.gm.chrgs.nsections);
#pragma omp critical
	    {
	      err = 0;
	    }	
	    //break;
	  }
	  pdens = ndens;
	}
      }// end master omp
    }// end parallel omp
  }// end else
  return err;
}
// -------------------- evolve openmp --------------------

// // -------------------- evolve openmp --------------------

// int NEvo::evolve_omp() {

//   int err;

//   // store the original nano delta time
//   double orig_dtn = nm.tm.ndeltat;
  
//   if(nm.rs.wsih4==1) {
//     BOOST_LOG_SEV(lg, info) << "With SiH4 rates";
//     err = compute_precompute_sih4();

//     // lsoda
//     ls_qsys = new ls_qsystem(this);

//     // IMPORTANT set nested parallelism
//     omp_set_nested(true);
  
//     clock_t begin_sim = std::clock();
//     ctime = 0.0;
//     err = -1;
//     int err_onestep = 0;
// #pragma omp parallel //shared(ctime, nm, err)
//     {//begin parallel omp
// #pragma omp master
//       {// begin master omp, controls the iterations in t (time)

// 	BOOST_LOG_SEV(lg, info) << "Set number of threads = " << omp_get_num_threads();
// 	BOOST_LOG_SEV(lg, info) << "Nested = " << omp_get_nested();
// 	BOOST_LOG_SEV(lg, info) << "Dynamic = " << omp_get_dynamic();
// 	BOOST_LOG_SEV(lg, info) << "Max active levels " << omp_get_max_active_levels();
// 	unsigned int npcount = 0;
// 	// get num_threads
// 	int num_threads = omp_get_num_threads();
// 	// execute loop in master thread
// 	double writing_step_begin = 0.0;

// 	// WARNING reset dtn
// 	nm.tm.ndeltat = orig_dtn;
// 	BOOST_LOG_SEV(lg, info) << "Reset to original dtn = " << nm.tm.ndeltat;
	
// 	for(unsigned int t=0; err<0; ++t, ++npcount) {
// 	  // update time
// 	  ctime += nm.tm.ndeltat;

// 	  //clock_t one_step_begin = std::clock();
// 	  double one_step_begin = omp_get_wtime();

// 	  //BOOST_LOG_SEV(lg, info) << "Number threads " << omp_get_num_threads();
// #pragma omp parallel /*num_threads(num_threads) shared(ctime)*/
// 	  {
// 	    err_onestep = evolve_one_step_omp(ctime);
// 	  }
// 	  if(npcount>static_cast<unsigned int>((0.01*nm.tm.tstop)/nm.tm.ndeltat)) {
// 	    clock_t elapsed_sim = std::clock();
// 	    double elapsed_secs = double(elapsed_sim - begin_sim) / CLOCKS_PER_SEC;
// 	    BOOST_LOG_SEV(lg, info) << "Writing results at simulation time = "
// 				    << ctime << " [elapsed secs = " 
// 				    << elapsed_secs << "]";
// 	    int err_ = write_partial_results(ctime);
// 	    //clock_t one_step_end = std::clock();
// 	    double one_step_end = omp_get_wtime();
// 	    //double elapsed_secs_onestep  = double(one_step_end - one_step_begin) / CLOCKS_PER_SEC;
// 	    double elapsed_secs_onestep = one_step_end - one_step_begin;
// 	    BOOST_LOG_SEV(lg, info) << "One step elapsed secs = "  << elapsed_secs_onestep;

// 	    double writing_step_end = omp_get_wtime();
// 	    double elapsed_secs_writing = writing_step_end - writing_step_begin;
// 	    BOOST_LOG_SEV(lg, info) << "Writing time/Delta elapsed sec = " << ctime 
// 				    << "\t" << elapsed_secs_writing;    
// 	    // reset clock
// 	    writing_step_begin = omp_get_wtime();
// 	    npcount=0;
// 	  }
// 	  if(ctime>nm.tm.tstop-nm.tm.ndeltat) {
// 	    BOOST_LOG_SEV(lg, info) << "Done!, time = " << ctime;
// 	    clock_t end_sim = std::clock();
// 	    double elapsed_secs = double(end_sim - begin_sim) / CLOCKS_PER_SEC;
// 	    BOOST_LOG_SEV(lg, info) << "Nanoparticle module elapsed time "
// 				    << elapsed_secs;
// 	    int err_ = write_partial_results(ctime);
// 	    err_ = write_dataset2d_hdf5(h5obj_nano, "Density", "density", ndens,
// 					cr.gm.vols.nsections, cr.gm.chrgs.nsections);
// #pragma omp critical
// 	    {
// 	      err = 0;
// 	    }	
// 	    //break;
// 	  }
// 	  pdens = ndens;
// 	}
//       }// end master omp
//     }// end parallel omp    
//   }
//   else {
//     // check = 2;
//     // // initial conditions, in boost vector_type x( 2 , 1.0 );
//     // N_Vector xini = NULL;
//     // xini = N_VNew_Serial(cr.gm.chrgs.nsections);
//     // NV_Ith_S(xini, 0) = RCONST(2.0);
//     // NV_Ith_S(xini, 1) = RCONST(1.0);

//     // std::cerr << "\n Initial condition : "
//     //           << NV_Ith_S(xini,0) << '\t' << NV_Ith_S(xini,1) << '\n' << NV_LENGTH_S(xini) << "\n"; 

//     // // Instantiate the solver
//     // sol = new Solver(this, xini, 0.0);
//     // //  Solver sol(xini, 0.0);
//     // check = 3;
//     // sol->print();
//     // // final time
//     // realtype tf = RCONST(50.0);
//     // // delta time
//     // realtype dt = RCONST(0.01);
//     // //   sol->compute(tf, dt);
//     // //   sol.compute(tf, dt);


//     // // or solve one step
//     // //   sol.compute_step(0.0, dt);

//     // // free memory for xini
//     // N_VDestroy_Serial(xini);

//     err = compute_precompute();

//     qsys = new qsystem(this);

//     auto stepper = odeint::make_controlled( 1.0e-6, 1.0e-6,
// 					    odeint::runge_kutta_dopri5< state_type >() );

  
//     // lsoda
//     ls_qsys = new ls_qsystem(this);

//     // IMPORTANT set nested parallelism
//     omp_set_nested(true);
  
//     clock_t begin_sim = std::clock();
//     ctime = 0.0;
//     err = -1;
// #pragma omp parallel //shared(ctime, nm, err)
//     {//begin parallel omp
// #pragma omp master
//       {// begin master omp, controls the iterations in t (time)

// 	BOOST_LOG_SEV(lg, info) << "Set number of threads = " << omp_get_num_threads();
// 	BOOST_LOG_SEV(lg, info) << "Nested = " << omp_get_nested();
// 	BOOST_LOG_SEV(lg, info) << "Dynamic = " << omp_get_dynamic();
// 	BOOST_LOG_SEV(lg, info) << "Max active levels " << omp_get_max_active_levels();
// 	unsigned int npcount = 0;
// 	// get num_threads
// 	int num_threads = omp_get_num_threads();
// 	// execute loop in master thread
// 	double writing_step_begin = 0.0;
// 	for(unsigned int t=0; err<0; ++t, ++npcount) {
// 	  // update time
// 	  ctime += nm.tm.ndeltat;

// 	  //clock_t one_step_begin = std::clock();
// 	  double one_step_begin = omp_get_wtime();

// 	  //BOOST_LOG_SEV(lg, info) << "Number threads " << omp_get_num_threads();
// #pragma omp parallel /*num_threads(num_threads) shared(ctime)*/
// 	  {
// 	    evolve_one_step_omp(ctime);
// 	  }
// 	  if(npcount>static_cast<unsigned int>((0.01*nm.tm.tstop)/nm.tm.ndeltat)) {
// 	    clock_t elapsed_sim = std::clock();
// 	    double elapsed_secs = double(elapsed_sim - begin_sim) / CLOCKS_PER_SEC;
// 	    BOOST_LOG_SEV(lg, info) << "Writing results at simulation time = "
// 				    << ctime << " [elapsed secs = " 
// 				    << elapsed_secs << "]";
// 	    int err_ = write_partial_results(ctime);
// 	    //clock_t one_step_end = std::clock();
// 	    double one_step_end = omp_get_wtime();
// 	    //double elapsed_secs_onestep  = double(one_step_end - one_step_begin) / CLOCKS_PER_SEC;
// 	    double elapsed_secs_onestep = one_step_end - one_step_begin;
// 	    BOOST_LOG_SEV(lg, info) << "One step elapsed secs = "  << elapsed_secs_onestep;

// 	    double writing_step_end = omp_get_wtime();
// 	    double elapsed_secs_writing = writing_step_end - writing_step_begin;
// 	    BOOST_LOG_SEV(lg, info) << "Writing time/Delta elapsed sec = " << ctime 
// 				    << "\t" << elapsed_secs_writing;    
// 	    // reset clock
// 	    writing_step_begin = omp_get_wtime();
// 	    npcount=0;
// 	  }
// 	  if(ctime>nm.tm.tstop-nm.tm.ndeltat) {
// 	    BOOST_LOG_SEV(lg, info) << "Done!, time = " << ctime;
// 	    clock_t end_sim = std::clock();
// 	    double elapsed_secs = double(end_sim - begin_sim) / CLOCKS_PER_SEC;
// 	    BOOST_LOG_SEV(lg, info) << "Nanoparticle module elapsed time "
// 				    << elapsed_secs;
// 	    int err_ = write_partial_results(ctime);
// 	    err_ = write_dataset2d_hdf5(h5obj_nano, "Density", "density", ndens,
// 					cr.gm.vols.nsections, cr.gm.chrgs.nsections);
// #pragma omp critical
// 	    {
// 	      err = 0;
// 	    }	
// 	    //break;
// 	  }
// 	  pdens = ndens;
// 	}
//       }// end master omp
//     }// end parallel omp
//   }// end else
//   return err;
// }


// // -------------------- evolve openmp --------------------

int NEvo::compute_initialdensity() {
  // initial density
  idens.resize(cr.grid);

  if (nm.ds.peakpos>cr.gm.vols.nsections) {
    BOOST_LOG_SEV(lg, fatal)
      << "Section number for peak of distribution is greater than max section"
      << " number " << nm.ds.peakpos << " > " << cr.gm.vols.nsections;
  }
  switch(nm.ds.distribution) {
    // Step distribution
    case 1:
    {
      int halfwidth = nm.ds.width / 2;
      for(unsigned int i=0; i<nm.ds.width; ++i) {
        idens[nm.ds.peakpos-halfwidth+i][cr.gm.chrgs.maxnegative] =
                              nm.ds.indens/static_cast<double>(nm.ds.width);
      }
      break;
    }
    // Gaussian distribution
    case 2:
    {
      // FIXME compute number density in sections to get Ntotal=sum Ni
      int halfwidth = nm.ds.width / 2;
      double vpeak  = cr.gm.vols.volumes[nm.ds.peakpos];
      double vsigma = vpeak
                    - cr.gm.vols.volumes[nm.ds.peakpos-halfwidth];
      for(unsigned int i=0; i<nm.ds.width; ++i) {
        unsigned int ix = static_cast<unsigned int>(nm.ds.peakpos-halfwidth+i);
        double v = cr.gm.vols.volumes[ix];

        idens[ix][cr.gm.chrgs.maxnegative] = gaussian_distribution(v, vpeak, vsigma);
      }
      break;
    }
    case 0:
//     default:
      idens[nm.ds.peakpos][cr.gm.chrgs.maxnegative] = nm.ds.indens;
      break;
  }

  if (nm.ds.chargewidth) {
    int pos = find_index(cr.gm.chrgs.charges, nm.ds.chargenegwidth);
    unsigned int negpos = 0;
    if (pos > 0) {
      negpos = static_cast<unsigned int>(pos);
    }
    else {
      negpos = 0;
      BOOST_LOG_SEV(lg, info) << "Warning not valid negative charge for width: "
			      << nm.ds.chargenegwidth;
    }
    BOOST_LOG_SEV(lg, info) << "Index, value for charge distribution width = "
			    << negpos << ", " << cr.gm.chrgs.charges[negpos];
    pos = find_index(cr.gm.chrgs.charges, nm.ds.chargeposwidth);

    unsigned int pospos = 0;
    if (pos > 0) {
      pospos = static_cast<unsigned int>(pos);
    }
    else {
      pospos = cr.gm.chrgs.nsections-1;
      BOOST_LOG_SEV(lg, info) << "Warning not valid positive charge for width: "
			      << nm.ds.chargeposwidth;
    }
    BOOST_LOG_SEV(lg, info) << "Index, value for charge distribution width = "
			    << pospos << ", " << cr.gm.chrgs.charges[pospos];

    double qwidth = static_cast<double>(cr.gm.chrgs.charges[pospos]
					-cr.gm.chrgs.charges[negpos]+1);

    darray i0dens(cr.gm.vols.nsections);
    for(unsigned int i=0; i<cr.gm.vols.nsections; ++i) {
      i0dens[i] = idens[i][cr.gm.chrgs.maxnegative]/qwidth;
    }
    
    for(unsigned int qpos = negpos; qpos < pospos+1; ++qpos) {
      for(unsigned int i=0; i<cr.gm.vols.nsections; ++i) {
        idens[i][qpos] = i0dens[i];
      }
    }
  }
  return 0;
}

int NEvo::compute_precompute() {

  int err;

  // Nanoparticle potential
  phid.resize(cr.grid);

  // Collision frequencies
  efreq.resize(cr.grid);
  ifreq.resize(cr.grid);
  tfreq.resize(cr.grid);

  compute_nanoparticle_potential();
  BOOST_LOG_SEV(lg, info) << "Writing nanoparticle potential to file.";
  err = write_dataset2d_hdf5(h5obj_nano, "Potential", "phid", phid,
                             cr.gm.vols.nsections, cr.gm.chrgs.nsections);

  compute_collisionfreq();

  err = write_collisionfreq();

  // initial density
  compute_initialdensity();
  
  // zero array
  zero2d.resize(cr.grid);

  for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
    for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
      zero2d[l][q] = 0.0;
    }
  }
  
  // auxiliar densities
  pdens.resize(cr.grid);
  pdens = idens;
  ndens.resize(cr.grid);

  pdens_aux.resize(cr.grid);
  ndens_aux.resize(cr.grid);

  // total rate
  trate2d.resize(cr.grid);
  dt2d.resize(cr.grid);
    
  // rate of coagulation surface growth charging nucleation
  crate2d.resize(cr.grid);
  srate2d.resize(cr.grid);
  qrate2d.resize(cr.grid);
  nrate2d.resize(cr.grid);

  srate2d_aux.resize(cr.grid);
  nrate2d_aux.resize(cr.grid);

  moments_file = new std::fstream(dirname + nano_filename + "-moments.dat", std::fstream::out);
  *moments_file << "#Time\tNumber\tVolume\tCharge\tmu11\tmu20\tmu02\tmu21\tmu12\tmu30\tmu03";

  moments.resize(10);
  compute_moments(idens);
  write_moments(0.0);

  // surface growth
  surface_area.resize(cr.gm.vols.nsections);
  surface_area = cr.gm.vols.radii.apply(area_from_radius<double>);
  surface_rate.resize(cr.gm.vols.nsections);
  surface_rate = nm.rs.sgrowth_rate * surface_area;

  adim_srate.resize(cr.gm.vols.nsections);

  adim_srate[0] = surface_rate[0] / (cr.gm.vols.volumes[1]-cr.gm.vols.volumes[0]);
  for (unsigned int l = 1; l < cr.gm.vols.nsections; ++l) {
    adim_srate[l] = surface_rate[l] / (cr.gm.vols.volumes[l+1]-cr.gm.vols.volumes[l]);
  }
  // Growth to last section only adds particles, removal requires to
  // add sections (note sign plus + )
  adim_srate[cr.gm.vols.nsections-1] = surface_rate[cr.gm.vols.nsections-2]
                           / (cr.gm.vols.volumes[cr.gm.vols.nsections-1]
                           -cr.gm.vols.volumes[cr.gm.vols.nsections-2]);

  // Write surface rate
  // TODO write m3/s growth rate surface_rate
  err = create_dataset_hdf5<darray>(h5obj_nano, adim_srate, "Rates",
                                    "surface_rate");
  err = create_attrib_hdf5(h5obj_nano, "/Rates", "srate_max", adim_srate.max());
  err = create_attrib_hdf5(h5obj_nano, "/Rates", "srate_min", adim_srate.min());

  gsurfacegrowth.resize(cr.grid);
  kcoagulation.resize(cr.grid);
  cfrequency.resize(cr.grid);
  jnucleation.resize(cr.grid);

  // vector of birth of particles in section
  birth_vector.resize(cr.grid);
  // vector of death of particles in section
  death_vector.resize(cr.grid);

  for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
    for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
      birth_vector[l][q]=0.0; death_vector[l][q]=0.0;
    }
  }
//   nqdens.resize(cr.gm.chrgs.nsections);

  return err;
}

// ============================== SiH4 ==============================

int NEvo::compute_precompute_sih4() {

  int err;

  // Nanoparticle potential
  phid.resize(cr.grid);

  // Collision frequencies
  efreq.resize(cr.grid);
  ifreq.resize(cr.grid);
  tfreq.resize(cr.grid);

  compute_nanoparticle_potential();
  BOOST_LOG_SEV(lg, info) << "Writing nanoparticle potential to file.";
  err = write_dataset2d_hdf5(h5obj_nano, "Potential", "phid", phid,
                             cr.gm.vols.nsections, cr.gm.chrgs.nsections);

  compute_collisionfreq();

  err = write_collisionfreq();

  // initial density
  compute_initialdensity();

  // zero array
  zero2d.resize(cr.grid);

  for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
    for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
      zero2d[l][q] = 0.0;
    }
  }
  
  // auxiliar densities
  pdens.resize(cr.grid);
  pdens = idens;
  ndens.resize(cr.grid);

  pdens_aux.resize(cr.grid);
  ndens_aux.resize(cr.grid);

  // total rate
  trate2d.resize(cr.grid);
  dt2d.resize(cr.grid);
  
  // rate of coagulation surface growth charging nucleation
  crate2d.resize(cr.grid);
  srate2d.resize(cr.grid);
  qrate2d.resize(cr.grid);
  nrate2d.resize(cr.grid);

  srate2d_aux.resize(cr.grid);
  nrate2d_aux.resize(cr.grid);

  moments_file = new std::fstream(dirname + nano_filename + "-moments.dat", std::fstream::out);
  *moments_file << "#Time\tNumber\tVolume\tCharge\tmu11\tmu20\tmu02\tmu21\tmu12\tmu30\tmu03";

  moments.resize(10);
  compute_moments(idens);
  write_moments(0.0);

  sgrowth_rate_sih4 = nm.rs.sgrowth_rate;
  // ==================== SiH4 ====================
  if(nm.rs.wsih4==1) {
    BOOST_LOG_SEV(lg, info) << "SiH4 parameters: ";
    sih4_file = new std::fstream(dirname + nano_filename + "-sih4.dat", std::fstream::out);
    *sih4_file << "#Time\tDensity\tSiH4Rate\tNucRate\tSGrowthRate";

    // SiH4 ratio x gas desity
    nsih4_ini = nm.rs.sih4ratio * pm.pars.neutral_density;
    nsih4 = nsih4_ini;
    nsih4_aux = nsih4;
    BOOST_LOG_SEV(lg, info) << "--> SiH4 initial density: " << nsih4;

    BOOST_LOG_SEV(lg, info) << "--> SiH4 mass density: " <<  cr.gm.gsys.nmdensity;

    // 1/rho_Si x mass SiH4
    sih4_vol = (1.0/cr.gm.gsys.nmdensity)*1.67e-27*28;
    BOOST_LOG_SEV(lg, info) << "--> SiH4 volume: " <<  sih4_vol;

    // SiH4 thermal velocity
    vsih4 = sqrt((8./pi)*(Kboltz*pm.pars.temperature)/nm.rs.sih4mass);
    BOOST_LOG_SEV(lg, info) << "--> SiH4 thermal velocity: " <<  vsih4;

    sgrowth_effective = nm.rs.sgrowth_rate/(nsih4_ini*vsih4*sih4_vol);
    BOOST_LOG_SEV(lg, info)
      << "--> SiH4 effective surface growth coefficient : "
      <<  sgrowth_effective;

    sgrowth_total_rate = 0.0;
  } 
  // surface growth
  surface_area.resize(cr.gm.vols.nsections);
  surface_area = cr.gm.vols.radii.apply(area_from_radius<double>);
  surface_rate.resize(cr.gm.vols.nsections);
  surface_rate = sgrowth_rate_sih4 * surface_area;

  adim_srate.resize(cr.gm.vols.nsections);

  adim_srate[0] = surface_rate[0] / (cr.gm.vols.volumes[1]-cr.gm.vols.volumes[0]);
  for (unsigned int l = 1; l < cr.gm.vols.nsections-1; ++l) {
    adim_srate[l] = surface_rate[l] / (cr.gm.vols.volumes[l+1]-cr.gm.vols.volumes[l]);
  }
  // Growth to last section only adds particles, removal requires to
  // add sections (note sign plus + )
  adim_srate[cr.gm.vols.nsections-1] = surface_rate[cr.gm.vols.nsections-2]
                           / (cr.gm.vols.volumes[cr.gm.vols.nsections-1]
                           -cr.gm.vols.volumes[cr.gm.vols.nsections-2]);

  // Write surface rate
  // TODO write m3/s growth rate surface_rate
  err = create_dataset_hdf5<darray>(h5obj_nano, adim_srate, "Rates",
                                    "surface_rate");
  err = create_attrib_hdf5(h5obj_nano, "/Rates", "srate_max", adim_srate.max());
  err = create_attrib_hdf5(h5obj_nano, "/Rates", "srate_min", adim_srate.min());

  gsurfacegrowth.resize(cr.grid);
  kcoagulation.resize(cr.grid);
  cfrequency.resize(cr.grid);
  jnucleation.resize(cr.grid);

  // vector of birth of particles in section
  birth_vector.resize(cr.grid);
  // vector of death of particles in section
  death_vector.resize(cr.grid);

  for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
    for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
      birth_vector[l][q]=0.0; death_vector[l][q]=0.0;
    }
  }
//   nqdens.resize(cr.gm.chrgs.nsections);

  return err;
}





int NEvo::compute_nanoparticle_potential() {
  for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
    for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
      phid[l][q] = Kcoul*cr.gm.chrgs.charges[q]*eCharge
                 *(1.0/cr.gm.vols.radii[l]);
    }
  }
  return 0;
}

int NEvo::compute_collisionfreq() {
  double max_efreq = 0.0;
  double max_ifreq = 0.0;
  double min_efreq = 1.0e15;
  double min_ifreq = 1.0e15;

  // From 1.Allen, J. E. Probe theory - the orbital motion approach. Phys. Scr. 45, 497 (1992).

  // FIXME optimization use a vector of areas instead of M_PI radii2[l]
  double kte = (2.0/3.0)*pm.es.emean*eCharge;
  //double kte = 3.0*pm.es.emean*eCharge;
  efreqfactor = 4.0 * M_PI * pm.es.ne * eCharge * sqrt(kte/(2.0*M_PI*eMass));

  BOOST_LOG_SEV(lg, info) << "Electron frequency prefactor: "
                          << efreqfactor;

  double ion_energy_from_temperature = (3.0/2.0) * Kboltz * pm.is.itemp;
  // WARNING ion_velocity not defined TODO
  double ion_velocity = 0.0;
  double ion_energy = ion_energy_from_temperature 
                    + 0.5*pm.is.imass*ion_velocity*ion_velocity;
  double kti = (2.0/3.0)*ion_energy;

  // ifreq factor
  ifreqfactor = 4.0 * M_PI * pm.is.ni * eCharge * sqrt(kti/(2.0*M_PI*pm.is.imass));

  BOOST_LOG_SEV(lg, info) << "Ion frequency prefactor: "
                          << ifreqfactor;

  darray radii2 = cr.gm.chrgs.charges*cr.gm.chrgs.charges;

//   std::cerr << "\n[ii] ifreqfactor" << ifreqfactor;
// #pragma omp parallel for collapse(2)
  for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
    for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
      efreq[l][q] = 0.0;
      ifreq[l][q] = 0.0;
      // negative Nanoparticles
      if (q<cr.gm.chrgs.maxnegative) {
        efreq[l][q] = efreqfactor * radii2[l]
                    * exp(eCharge*phid[l][q]/kte);
        ifreq[l][q] = ifreqfactor * radii2[l]
                    * (1.0 - eCharge*phid[l][q]/kti);
      } else {
        // positive Nanoparticles
        if (q>cr.gm.chrgs.maxnegative) {
          efreq[l][q] = efreqfactor * radii2[l]
                      * (1.0 + eCharge*phid[l][q]/kte);
          ifreq[l][q] = ifreqfactor * radii2[l]
                      * exp(-eCharge*phid[l][q]/kti);
        } else { // Neutral
          efreq[l][q] = efreqfactor * radii2[l];
          ifreq[l][q] = ifreqfactor * radii2[l];
        }
      }
      if (efreq[l][q] > max_efreq) {
        max_efreq = efreq[l][q];
      }
      if (ifreq[l][q] > max_ifreq) {
        max_ifreq = ifreq[l][q];
      }
      if (efreq[l][q] < min_efreq) {
        min_efreq = efreq[l][q];
      }
      if (ifreq[l][q] < min_ifreq) {
        min_ifreq = ifreq[l][q];
      }
    }

  }

  double max_tfreq = 0.0;
  double min_tfreq = 1.0e20;

  if(nm.nano.tunnel == 1) {
    for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
      for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
        tfreq[l][q] = 0.0;
        // WARNING tunnel current only possible when electrons are in the particle
        if(q<cr.gm.chrgs.maxnegative){
          tfreq[l][q] = ptunnel(cr.gm.chrgs.charges[q], cr.gm.vols.radii[l]);
          if (!std::isfinite(tfreq[l][q])) {
            tfreq[l][q] = 0.0;
          }
//           std::cerr << std::endl << "[ii] tfreq " << tfreq[l][q] << std::endl;
  //         if ((tfreq[l][q] > max_efreq)/*&&(tfreq[l][q] > max_ifreq)*/) {
  //            tfreq[l][q] = max_efreq;
  //         }
  //         if (tfreq[l][q] > max_tfreq) {
  //           max_tfreq = tfreq[l][q];
  //         }
  //         if (tfreq[l][q] < min_tfreq) {
  //           min_tfreq = tfreq[l][q];
  //         }
  //         if ((tfreq[l][q] > max_efreq)/*&&(tfreq[l][q] > max_ifreq)*/) {
  //            tfreq[l][q] = max_efreq;
  //         }
          if (tfreq[l][q] > max_tfreq) {
            max_tfreq = tfreq[l][q];
          }
          if (tfreq[l][q] < min_tfreq) {
            min_tfreq = tfreq[l][q];
          }
//           if ((tfreq[l][q] > efreq[l][q])/*&&(tfreq[l][q] > max_ifreq)*/) {
//             tfreq[l][q] = efreq[l][q];
//           }
          //WARNING if tunnel frequency is too higher
          // FIXME these allow to accelerate the calculations
          if ((tfreq[l][q] > efreq[l][q])/*&&(tfreq[l][q] > max_ifreq)*/) {
             efreq[l][q] = 0.0;
             tfreq[l][q] = 0.0;
          }
        }
      }
    }
  }


  BOOST_LOG_SEV(lg, info) << "Max electron collision frequency "
            << max_efreq;
  BOOST_LOG_SEV(lg, info) << "Max ion collision frequency "
            << max_ifreq;
  BOOST_LOG_SEV(lg, info) << "Max electron + ion collision frequency "
            << max_efreq + max_ifreq;
  BOOST_LOG_SEV(lg, info) << "Max tunnel collision frequency "
            << max_tfreq;
  BOOST_LOG_SEV(lg, info) << "Min electron collision frequency "
            << min_efreq;
  BOOST_LOG_SEV(lg, info) << "Min ion collision frequency "
            << min_ifreq;
  BOOST_LOG_SEV(lg, info) << "Min electron + ion collision frequency "
            << min_efreq + min_ifreq;
  BOOST_LOG_SEV(lg, info) << "Min tunnel collision frequency "
            << min_tfreq;

  return 0;
}

int NEvo::write_collisionfreq() {
  int err;

  BOOST_LOG_SEV(lg, info) << "Writing electron frequency to file.";
  err = write_dataset2d_hdf5(h5obj_nano, "Collision_frequencies",
                             "efreq",
                             efreq, cr.gm.vols.nsections,
                             cr.gm.chrgs.nsections);
  // Write charging rate
  double rele_max = *std::max_element(efreq.origin(), efreq.origin()
                                                      + efreq.num_elements());
  err = create_attrib_hdf5(h5obj_nano, "/Collision_frequencies", "rele_max", rele_max);
  double rele_min = *std::min_element(efreq.origin(), efreq.origin()
                                                      + efreq.num_elements());
  err = create_attrib_hdf5(h5obj_nano, "/Collision_frequencies", "rele_min", rele_min);  

  BOOST_LOG_SEV(lg, info) << "Writing ion frequency to file.";
  err = write_dataset2d_hdf5(h5obj_nano, "Collision_frequencies",
                             "ifreq",
                             ifreq, cr.gm.vols.nsections,
                             cr.gm.chrgs.nsections);
  // Write charging rate
  double rion_max = *std::max_element(ifreq.origin(), ifreq.origin()
                                                      + ifreq.num_elements());
  err = create_attrib_hdf5(h5obj_nano, "/Collision_frequencies", "rion_max", rion_max);
  double rion_min = *std::min_element(ifreq.origin(), ifreq.origin()
                                                      + ifreq.num_elements());
  err = create_attrib_hdf5(h5obj_nano, "/Collision_frequencies", "rion_min", rion_min);


  BOOST_LOG_SEV(lg, info) << "Writing tunnel frequency to file.";
  err = write_dataset2d_hdf5(h5obj_nano, "Collision_frequencies",
                             "tfreq",
                             tfreq, cr.gm.vols.nsections,
                             cr.gm.chrgs.nsections);
  // Write charging rate
  double rtun_max = *std::max_element(tfreq.origin(), tfreq.origin()
                                                      + tfreq.num_elements());
  err = create_attrib_hdf5(h5obj_nano, "/Collision_frequencies", "rtun_max", rtun_max);
  double rtun_min = *std::min_element(tfreq.origin(), tfreq.origin()
                                                      + tfreq.num_elements());
  err = create_attrib_hdf5(h5obj_nano, "/Collision_frequencies", "rtun_min", rtun_min);

  return err;
}

int NEvo::compute_moments(boost_array2d dens, darray &moments) {
  darray vpivots = cr.gm.vols.volumes;
  darray qpivots = cr.gm.chrgs.charges;
  moments = 0.0;
  for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
    for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
      moments[0] += dens[l][q];              //total_number
      moments[1] += dens[l][q] * vpivots[l]; //total_volume
      moments[2] += dens[l][q] * qpivots[q]; //total_charge
      // moments[3] += dens[l][q] * vpivots[l] * qpivots[q];//mu_11
      // moments[4] += dens[l][q] * vpivots[l] * vpivots[l];//mu_20
      // moments[5] += dens[l][q] * qpivots[q] * qpivots[q];//mu_02
      // moments[6] += dens[l][q] * vpivots[l] * vpivots[l] * qpivots[q];//mu_21
      // moments[7] += dens[l][q] * vpivots[l] * qpivots[q] * qpivots[q];//mu_12
      // moments[8] += dens[l][q] * vpivots[l] * vpivots[l] * vpivots[l];//mu_30
      // moments[9] += dens[l][q] * qpivots[q] * qpivots[q] * qpivots[q];//mu_03
    }
  }
  return 0;
}

int NEvo::compute_moments(boost_array2d dens) {
  darray vpivots = cr.gm.vols.volumes;
  darray qpivots = cr.gm.chrgs.charges;
  moments = 0.0;
  for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
    for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
      moments[0] += dens[l][q];              //total_number
      moments[1] += dens[l][q] * vpivots[l]; //total_volume
      moments[2] += dens[l][q] * qpivots[q]; //total_charge
      moments[3] += dens[l][q] * vpivots[l] * qpivots[q];//mu_11
      moments[4] += dens[l][q] * vpivots[l] * vpivots[l];//mu_20
      moments[5] += dens[l][q] * qpivots[q] * qpivots[q];//mu_02
      moments[6] += dens[l][q] * vpivots[l] * vpivots[l] * qpivots[q];//mu_21
      moments[7] += dens[l][q] * vpivots[l] * qpivots[q] * qpivots[q];//mu_12
      moments[8] += dens[l][q] * vpivots[l] * vpivots[l] * vpivots[l];//mu_30
      moments[9] += dens[l][q] * qpivots[q] * qpivots[q] * qpivots[q];//mu_03
    }
  }
  return 0;
}

int NEvo::write_moments(double ctime) {
  *moments_file << '\n' << ctime;
  for(unsigned int i=0; i<moments.size(); ++i) {
    *moments_file << '\t' << moments[i];
  }
  moments_file->flush();
//   if(moments[0]!=moments[0]) std::terminate;
  if(!std::isfinite(moments[0])){
    std::cerr << "\n[ee] Particle number is not finite. Terminate.\n";
    BOOST_LOG_SEV(lg, fatal) << "\nParticle number is not finite. Terminate.\n";
    std::terminate();
  }
  if(sih4_file != NULL) {
    *sih4_file << '\n' << ctime
	       << '\t' << nsih4
	       << '\t' << sih4rate
               << '\t' << nm.rs.nucleation_rate*nsih4/nsih4_ini
	       << '\t' << sgrowth_rate_sih4
	       << '\t' << sgrowth_total_rate;
    sih4_file->flush();
  }
  return 0;
}

int NEvo::compute_explicit_charging(double dt) {

// dynamic plasma must recompute collisionfreq
//   compute_collisionfreq();

  double tun = 0.0;
  if(nm.nano.tunnel==1) {
    tun = 1.0;
  }
// #pragma omp for collapse(2)
  for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
    for (unsigned int q = 1; q < cr.gm.chrgs.nsections-1; ++q) {

//       ndens[l][q] = pdens[l][q] *(1.0 - dt*((ifreq[l][q]+tun*tfreq[l][q]) + efreq[l][q]))
//                   + dt *((ifreq[l][q-1]+tun*tfreq[l][q-1])*pdens[l][q-1]
//                         + efreq[l][q+1]*pdens[l][q+1]);
      double hdt = dt;
      double dens;
      do {
        for(double tp=0.0; tp<dt; tp+=hdt){
          dens = pdens[l][q] *(1.0 - hdt*((ifreq[l][q]+tun*tfreq[l][q]) + efreq[l][q]))
                + hdt *((ifreq[l][q-1]+tun*tfreq[l][q-1])*pdens[l][q-1]
                      + efreq[l][q+1]*pdens[l][q+1]);
        }
        hdt *= 0.25;
      } while(dens<0.0);
      ndens[l][q] = dens;
    }
  }

  for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
//     ndens[l][0] = pdens[l][0]
//                          + dt *(efreq[l][1]*pdens[l][1]
//                               - (ifreq[l][0]+tun*tfreq[l][0])*pdens[l][0]);

    double hdt = dt;
    double dens;
    do {
      for(double tp=0.0; tp<dt; tp+=hdt){
        dens = pdens[l][0]
             + hdt *(efreq[l][1]*pdens[l][1]
                    - (ifreq[l][0]+tun*tfreq[l][0])*pdens[l][0]);
      }
      hdt *= 0.25;
    } while(dens<0.0);
    ndens[l][0] = dens;

    hdt = dt;

    do {
      for(double tp=0.0; tp<dt; tp+=hdt) {
        dens = pdens[l][cr.gm.chrgs.nsections-1]
                      + hdt *((ifreq[l][cr.gm.chrgs.nsections-2]
                            +tun*tfreq[l][cr.gm.chrgs.nsections-2])
                              *pdens[l][cr.gm.chrgs.nsections-2]
                            - efreq[l][cr.gm.chrgs.nsections-1]
                            *pdens[l][cr.gm.chrgs.nsections-1]);
      }
      hdt *= 0.5;
    } while(dens<0.0);
    ndens[l][cr.gm.chrgs.nsections-1] = dens;

  }

//   bool negative_density = false;
// // #pragma omp for collapse(2)
//   for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
//     for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
// 
//       if (ndens[l][q]<0.0){//EPSILON) {
// //         std::cout << "\n[ee] Negative nanoparticle charging =\t"
// //                   << ndens[l][q] << " " << l << " " << q << " .Abort.\n";
// //         std::cout << "radii = " << radii[l] << "\n charge = " << chrgs.pivots[q] << "\n";
// //         BOOST_LOG_SEV(lg, fatal) << "Negative nanoparticle charging, ndens[" << l << "]["
// //                   << q << "] = " << ndens[l][q];
// //         BOOST_LOG_SEV(lg, info) << "radius = " << cr.gm.vols.radii[l]
// //                                 << "\t charge = " << cr.gm.chrgs.charges[q];
// //         BOOST_LOG_SEV(lg, info) << "ifreq = " << ifreq[l][q]
// //                                 << "\t efreq = " << efreq[l][q]
// //                                 << "tfreq = " << tfreq[l][q];
// // 
// //         std::cerr << "\n[ee] Negative nanoparticle charging, ndens[" << l << "]["
// //                   << q << "] = " << ndens[l][q] << '\n';
//            negative_density = true;       
// //         std::terminate();
//         ndens[l][q]=0.0;
//       }
//     }
//   }
//   
//   if(negative_density) {
//     BOOST_LOG_SEV(lg, fatal) << "Negative nanoparticle density, charging";
//   }
//   std::cout << "\n[ii] ntot ini = " << tot;
//   tot = 0.0;
//   for (unsigned int l = 0; l < vols.nsecs; ++l) {
//     for (unsigned int q = qind[l][0]; q < qind[l][1]+1; ++q) {
//       tot += ndens[l][q] * vpivots[l];
//     }
//   }
//   std::cout << "\n[ii] ntot fin = " << tot;
//   for (unsigned int l = 0; l < vols.nsecs; l++) {
//     for (unsigned int q = qind[l][0]; q < qind[l][1]+1; q++) {
//       pdens[l][q] = ndens[l][q];
//     }
//   }
  return 0;
}


int NEvo::write_one_step(double ctime,
			 double begin_sim,
			 double one_step_begin) {

  clock_t elapsed_sim = std::clock();
  double elapsed_secs = double(elapsed_sim - begin_sim) / CLOCKS_PER_SEC;
  BOOST_LOG_SEV(lg, info) << "Writing results at simulation time = "
			  << ctime << " [elapsed secs = " 
			  << elapsed_secs << "]";
  int err_ = write_partial_results(ctime);
  double one_step_end = omp_get_wtime();
  double elapsed_secs_onestep = one_step_end - one_step_begin;
  BOOST_LOG_SEV(lg, info) << "One step elapsed secs = "  << elapsed_secs_onestep;

  return 0;
}

int NEvo::write_partial_results(double ctime) {

  std::string density_name("density"+std::to_string(ctime));
  
  int err;
  err = write_dataset2d_hdf5(h5obj_nano, "Density", density_name, ndens,
                             cr.gm.vols.nsections, cr.gm.chrgs.nsections);

  // Current time of simulation
  err = create_dsattrib_hdf5(h5obj_nano, "Density", density_name, "time", ctime);
  err = create_attrib_hdf5(h5obj_nano, "Time", density_name, ctime);

  // Charging rate
  std::string rate_name("rate"+std::to_string(ctime));
  err = write_dataset2d_hdf5(h5obj_nano, "Charging", rate_name, qrate2d,
                             cr.gm.vols.nsections, cr.gm.chrgs.nsections);
  err = create_dsattrib_hdf5(h5obj_nano, "Charging", rate_name, "time", ctime);
  
  // Coagulation rate
  err = write_dataset2d_hdf5(h5obj_nano, "Coagulation", rate_name, crate2d,
                             cr.gm.vols.nsections, cr.gm.chrgs.nsections);
  err = create_dsattrib_hdf5(h5obj_nano, "Coagulation", rate_name, "time", ctime);
  
  // Surface growth rate
  err = write_dataset2d_hdf5(h5obj_nano, "Surface_growth", rate_name, srate2d,
                             cr.gm.vols.nsections, cr.gm.chrgs.nsections);
  err = create_dsattrib_hdf5(h5obj_nano, "Surface_growth", rate_name, "time", ctime);
  
  // Nucleation rate
  err = write_dataset2d_hdf5(h5obj_nano, "Nucleation", rate_name, nrate2d,
                             cr.gm.vols.nsections, cr.gm.chrgs.nsections);
  err = create_dsattrib_hdf5(h5obj_nano, "Nucleation", rate_name, "time", ctime);
  
  compute_moments(ndens);
  write_moments(ctime);

  return 0;
}

int NEvo::evolve_one_step(double ctime) {
  int err;

  state_type nqdens(cr.gm.chrgs.nsections);
  
//   std::vector<double> nqdens(cr.gm.chrgs.nsections);


  double sdt=0.0;
  if(nm.rs.wch == 1) {
    clock_t begin_tcharge = std::clock();

//     for(;;) {
//       if(sdt>nm.tm.ndeltat) break;
//       err = compute_explicit_charging(nm.tm.qdeltat);
//   //     std::cerr << sdt;
//       sdt += nm.tm.qdeltat;
//       pdens = ndens;
//     }
    
    

    for(unsigned int l=0; l<cr.gm.vols.nsections; ++l) {
//       auto stepper = odeint::make_controlled( 1.0e-3, 1.0e-3,
//                                           odeint::runge_kutta_dopri5< state_type >() );
      auto stepper = odeint::make_controlled( 1.0e-6, 1.0e-6,
                                             odeint::runge_kutta_cash_karp54< state_type >() );
      qsys->l = l;
//       std::vector<double> nqdens(cr.gm.chrgs.nsections);
      for(unsigned int q=0; q<cr.gm.chrgs.nsections; ++q) {
         nqdens[q] = pdens[l][q];
      }
      
//       bool success = false;
      double ldtq = nm.tm.qdeltat;
      double ldtn = nm.tm.ndeltat;
      double ttime = ctime;

//       qsystem qs(qsys);
      odeint::integrate_adaptive(stepper_type(), *qsys, nqdens, ttime,
                                 ttime+ldtn, ldtq);
//       odeint::integrate(*qsys, nqdens, ttime,
//                                  ttime+ldtn, ldtq);

      for(unsigned int q=0; q<cr.gm.chrgs.nsections; ++q) {
        ndens[l][q] = nqdens[q];
      }
    }
    pdens = ndens;
    clock_t end_tcharge = std::clock();
    double elapsed_secs = double(end_tcharge - begin_tcharge) / CLOCKS_PER_SEC;
    BOOST_LOG_SEV(lg, info) << "Charging elapsed secs = " << elapsed_secs;
  }
  
  clock_t begin_tgrowth = std::clock();

  err = advance_nocharging(ctime);

  clock_t end_tgrowth = std::clock();
  double elapsed_secs = double(end_tgrowth - begin_tgrowth) / CLOCKS_PER_SEC;
  BOOST_LOG_SEV(lg, info) << "Growth elapsed secs = " << elapsed_secs;
  
  return err;
}


// ------------------ openmp version ------------------

// ==================== LSODA

int NEvo::evolve_one_step_omp(double ctime) {

  int err, reterr = 0;

  darray new_mom(0.0, 3);
  darray old_mom(0.0, 3);
  //mom.resize(3);

#pragma omp master
  {
    if(nm.rs.wch == 1) {
    INIT:
      // compute_moments(pdens, new_mom);
	
      ls_qsystem qsysy(*ls_qsys);
      double wtime = omp_get_wtime();
      double min_dtq = 1.0;
#pragma omp parallel firstprivate(qsysy)// shared(min_dtq, ctime, ndens, pdens, qrate2d)//, ctime, qsys)
      {
#pragma omp for nowait schedule(dynamic)// ordered schedule(static, 1)// nowait
	for(unsigned int l=0; l<cr.gm.vols.nsections; ++l) {
	  // work with section l
	  qsysy.l = l;
	  // each thread owns a density vector (fixed l)
	  state_type nqdens(cr.gm.chrgs.nsections);
	  // tolerances for lsoda
	  double atol[cr.gm.chrgs.nsections], rtol[cr.gm.chrgs.nsections];
	  // set the target density
	  for(unsigned int q=0; q<cr.gm.chrgs.nsections; ++q) {
	    nqdens[q] = pdens[l][q];
	    atol[q] = 1.0e-4;
	    rtol[q] = 1.0e-2;
	    //sum1 += pdens[l][q];
	  }
	  // get time
	  double ldtq = nm.tm.qdeltat;
	  double ldtn = nm.tm.ndeltat;
	  double ttime = ctime;	  
	  // INTEGRATE LSODA begins
	  // set options
	  struct lsoda_opt_t opt = {0};
	  //opt.h0 = ldtq;                // initial time step
	  opt.ixpr = 0;                 // additional printing
	  opt.rtol = rtol;              // relative tolerance
	  opt.atol = atol;              // absolute tolerance
	  opt.itask = 1;                // normal integration
	  opt.mxstep = 1000000000;      // max steps
	  // set lsoda context
	  struct lsoda_context_t ctx = {
					.function = ls_qsystemf,
					.data = &qsysy,
					.neq = static_cast<int>(cr.gm.chrgs.nsections),
					.state = 1,
	  };
	  lsoda_prepare(&ctx, &opt);
	  // time for lsoda
	  double tin = ttime;
	  double tout = ttime+ldtn;
	  // integrate
	  lsoda(&ctx, nqdens.data(), &tin, tout);
	  min_dtq = std::min(min_dtq, ctx.common->hu);
	  if (ctx.state <= 0) {
	    std::cerr << "\nerror istate = " << ctx.state;
	    //std::terminate();
	  }
	  // free context
	  lsoda_free(&ctx);
	  
	  for(unsigned int q=0; q<cr.gm.chrgs.nsections; ++q) {
	    ndens[l][q] = (nqdens[q]>0? nqdens[q]: 0.0);
	    qrate2d[l][q] = (ndens[l][q]-pdens[l][q])/ldtn;
	  }	  
	}//for l in sections
      }// parallel region
      //#pragma omp flush

      // compute moments after charging
      // compute_moments(ndens, old_mom);
      // BOOST_LOG_SEV(lg, info) << "number rel error " << ctime << '\t' << relative_error<double>(old_mom[0], new_mom[0]);
      // BOOST_LOG_SEV(lg, info) << "number abs error " << ctime << '\t' << absolute_error<double>(old_mom[0], new_mom[0]);
      // BOOST_LOG_SEV(lg, info) << "volume rel error " << ctime << '\t' << relative_error<double>(old_mom[1], new_mom[1]);
      // BOOST_LOG_SEV(lg, info) << "volume abs error " << ctime << '\t' << absolute_error<double>(old_mom[1], new_mom[1]);
     
      
      // if((is_close(new_mom[0], old_mom[0])==false) ||
      // 	 (is_close(new_mom[1], old_mom[1])==false)) {
      // 	BOOST_LOG_SEV(lg, info) << "Error in charging, moments differ";
      // 	BOOST_LOG_SEV(lg, info) << "pdens N = " << old_mom[0];
      // 	BOOST_LOG_SEV(lg, info) << "ndens N = " << new_mom[0];
      // 	BOOST_LOG_SEV(lg, info) << "cn0 " << abs(new_mom[0]-old_mom[0]);
      //   BOOST_LOG_SEV(lg, info) << "pdens V = " << old_mom[1];
      //   BOOST_LOG_SEV(lg, info) << "ndens V = " << new_mom[1];
      // 	BOOST_LOG_SEV(lg, info) << "mv0 " << abs(new_mom[1]-old_mom[1]);
      // 	//std::terminate();
      // }
      
      pdens = ndens;

      // wtime = omp_get_wtime() - wtime;
      // BOOST_LOG_SEV(lg, info) << "Charging elapsed secs = " << wtime;
      // BOOST_LOG_SEV(lg, info) << "Min dtq " << min_dtq;
    }//charging block    
    
    // growth clock
    double begin_tgrowth = omp_get_wtime();
    if(nm.rs.wsih4==1) {
      //err = advance_nocharging_ompadpsih4(ctime);
      err = advance_nosplit_ompadpsih4(ctime);
      if(err == -1) {
	reterr = -1;
	nm.tm.ndeltat *= 0.95;
	BOOST_LOG_SEV(lg, info) << "t, New dt = " << ctime << '\t' << nm.tm.ndeltat;
	goto INIT;
      }
      else {
	reterr = 0;
      }
    }
    else {
      err = advance_nocharging_ompadp(ctime);
    }
    double end_tgrowth = omp_get_wtime();
    double elapsed_secs = end_tgrowth - begin_tgrowth;
    BOOST_LOG_SEV(lg, info)
      << "Growth (coagulation+sgrowth+nucleation) elapsed secs = "
      << elapsed_secs;
    
  }//omp master
  return reterr;
}

// full adaptive lsoda

// ==================== LSODA

int NEvo::evolve_one_step_adapt(double ctime) {

  int err, reterr = 0;

  darray new_mom(0.0, 3);
  darray old_mom(0.0, 3);
  //mom.resize(3);

#pragma omp master
  {
    if(nm.rs.wch == 1) {
    INIT:
      // compute_moments(pdens, new_mom);
	
      ls_qsystem qsysy(*ls_qsys);
      double wtime = omp_get_wtime();
      double min_dtq = 1.0;
#pragma omp parallel firstprivate(qsysy)// shared(min_dtq, ctime, ndens, pdens, qrate2d)//, ctime, qsys)
      {
#pragma omp for nowait schedule(dynamic)// ordered schedule(static, 1)// nowait
	for(unsigned int l=0; l<cr.gm.vols.nsections; ++l) {
	  // work with section l
	  qsysy.l = l;
	  // each thread owns a density vector (fixed l)
	  state_type nqdens(cr.gm.chrgs.nsections);
	  // tolerances for lsoda
	  double atol[cr.gm.chrgs.nsections], rtol[cr.gm.chrgs.nsections];
	  // set the target density
	  for(unsigned int q=0; q<cr.gm.chrgs.nsections; ++q) {
	    nqdens[q] = pdens[l][q];
	    atol[q] = 1.0e-4;
	    rtol[q] = 1.0e-2;
	    //sum1 += pdens[l][q];
	  }
	  // get time
	  double ldtq = nm.tm.qdeltat;
	  double ldtn = nm.tm.ndeltat;
	  double ttime = ctime;	  
	  // INTEGRATE LSODA begins
	  // set options
	  struct lsoda_opt_t opt = {0};
	  //opt.h0 = ldtq;                // initial time step
	  opt.ixpr = 0;                 // additional printing
	  opt.rtol = rtol;              // relative tolerance
	  opt.atol = atol;              // absolute tolerance
	  opt.itask = 1;                // normal integration
	  opt.mxstep = 1000000000;      // max steps
	  // set lsoda context
	  struct lsoda_context_t ctx = {
					.function = ls_qsystemf,
					.data = &qsysy,
					.neq = static_cast<int>(cr.gm.chrgs.nsections),
					.state = 1,
	  };
	  lsoda_prepare(&ctx, &opt);
	  // time for lsoda
	  double tin = ttime;
	  double tout = ttime+ldtn;
	  // integrate
	  lsoda(&ctx, nqdens.data(), &tin, tout);
	  min_dtq = std::min(min_dtq, ctx.common->hu);
	  if (ctx.state <= 0) {
	    std::cerr << "\nerror istate = " << ctx.state;
	    //std::terminate();
	  }
	  // free context
	  lsoda_free(&ctx);
	  
	  for(unsigned int q=0; q<cr.gm.chrgs.nsections; ++q) {
	    ndens[l][q] = (nqdens[q]>0? nqdens[q]: 0.0);
	    qrate2d[l][q] = (ndens[l][q]-pdens[l][q])/ldtn;
	  }	  
	}//for l in sections
      }// parallel region
      //#pragma omp flush    
      pdens = ndens;

    }//charging block    
    
    // growth clock
    double begin_tgrowth = omp_get_wtime();
    if(nm.rs.wsih4==1) {
      //err = advance_nocharging_ompadpsih4(ctime);
      err = advance_nosplit_sih4_adapt(ctime);
      if(err == -1) {
	reterr = -1;
	nm.tm.ndeltat *= 0.95;
	BOOST_LOG_SEV(lg, info) << "t, New dt = " << ctime << '\t' << nm.tm.ndeltat;
	goto INIT;
      }
      else {
	reterr = 0;
      }
    }
    else {
      err = advance_nocharging_adapt(ctime);
    }
    double end_tgrowth = omp_get_wtime();
    double elapsed_secs = end_tgrowth - begin_tgrowth;
    // BOOST_LOG_SEV(lg, info)
    //   << "Growth (coagulation+sgrowth+nucleation) elapsed secs = "
    //   << elapsed_secs;    
    
  }//omp master
  return reterr;
}



// ++++++++++++++++++++ BOOST working version

// int NEvo::evolve_one_step_omp(double ctime) {

//   int err;

// #pragma omp master
//   {
//     //state_type nqdens(cr.gm.chrgs.nsections);
//     double wtime = omp_get_wtime();
//     double sdt=0.0;
//     if(nm.rs.wch == 1) {
//       //clock_t begin_tcharge = std::clock();
//       //BOOST_LOG_SEV(lg, info) << "Computing charges";
//       #pragma omp parallel default(none) shared(ctime, ndens, pdens, qsys, qrate2d)//, ctime, qsys)
//       {
// 	//std::cout << "\n[ii] number of threads in the team " << omp_get_thread_num();
// 	// dynamic scheduling because the load is unbalanced
// 	// the odeint solver is fast for density = 0 but takes time when
// 	// density != 0
// 	#pragma omp for nowait schedule(dynamic)// ordered schedule(static, 1)// nowait
// 	for(unsigned int l=0; l<cr.gm.vols.nsections; ++l) {
// 	  //       auto stepper = odeint::make_controlled( 1.0e-3, 1.0e-3,
// 	  //                                           odeint::runge_kutta_dopri5< state_type >() );
// 	  auto stepper = odeint::make_controlled( 1.0e-3, 1.0e-3,
// 						  odeint::runge_kutta_cash_karp54< state_type >() );
// 	  //auto stepper = odeint::make_controlled( 1.0e-3, 1.0e-3,
// 	  //					odeint::runge_kutta_fehlberg78< state_type >() );
// 	  //runge_kutta_fehlberg78
// 	  //qsys->l = l;

// 	  //wtime = omp_get_wtime();
// 	  // make a copy of qsys for each thread
// 	  qsystem qsysy(*qsys);
// 	  // work with section l
// 	  qsysy.l = l;
// 	  // each thread owns a density vector (fixed l)
// 	  state_type nqdens(cr.gm.chrgs.nsections);
// 	  // wtime = omp_get_wtime() - wtime;
// 	  // BOOST_LOG_SEV(lg, info) << "Copy elapsed secs = " << wtime;
// 	  // wtime = omp_get_wtime();

// 	  //double sum1 = 0;
// 	  for(unsigned int q=0; q<cr.gm.chrgs.nsections; ++q) {
// 	    nqdens[q] = pdens[l][q];
// 	    //sum1 += pdens[l][q];
// 	  }
// 	  //std::cout << "\n[ii] thread number " << omp_get_thread_num() << " l :  " << qsys->l << " sum1 " << sum1;
	  
// 	  //       bool success = false;
// 	  double ldtq = nm.tm.qdeltat;
// 	  double ldtn = nm.tm.ndeltat;
// 	  double ttime = ctime;
	  
// 	  odeint::integrate_adaptive(stepper_type(), qsysy, nqdens, ttime,
// 	  			     ttime+ldtn, ldtq);
// 	  // wtime = omp_get_wtime() - wtime;
// 	  // BOOST_LOG_SEV(lg, info) << "Solver elapsed secs = " << wtime;
// 	  //odeint::integrate_const( odeint::runge_kutta4< state_type >() ,qsysy, nqdens, ttime,
// 	  //			  ttime+ldtn, ldtq);
	  
// 	  //double sum = 0;
// 	  for(unsigned int q=0; q<cr.gm.chrgs.nsections; ++q) {
// 	    ndens[l][q] = nqdens[q];
// 	    //sum += nqdens[q];
// 	    //std::cout << "\n[ii] thread number " << omp_get_thread_num() << " l, q " << qsys->l << ", " << q;
// 	    qrate2d[l][q] = (ndens[l][q]-pdens[l][q])/ldtn;
// 	  }
// 	  //std::cout << "\n[ii] thread number " << omp_get_thread_num() << " l :  " << qsysy.l << " sum " << sum;
// 	  //std::cout << "\n[ii] thread number " << omp_get_thread_num() << " TIME " << ttime << "\t" << ctime;
// 	  //wtime = omp_get_wtime() - wtime;
// 	  //std::cout << "\n[ii] thread number " << omp_get_thread_num() << " elapsed time " << wtime;
//       }//for l in sections
//     }// parallel region
//     //#pragma omp flush
//     pdens = ndens;
//     //clock_t end_tcharge = std::clock();
//     //double elapsed_secs = double(end_tcharge - begin_tcharge) / CLOCKS_PER_SEC;
//     //BOOST_LOG_SEV(lg, info) << "Charging elapsed secs = " << elapsed_secs;
//   }//charging block
  
//   wtime = omp_get_wtime() - wtime;
//   BOOST_LOG_SEV(lg, info) << "Charging elapsed secs = " << wtime;

//   // growth clock
//   wtime = omp_get_wtime() - wtime;
//   clock_t begin_tgrowth = std::clock();
  
// //#pragma omp parallel num_threads(10)
// //    {
//     err = advance_nocharging_ompadp(ctime);
//     //BOOST_LOG_SEV(lg, info) << "No charging thread " << omp_get_thread_num();
//     //    }
    
//     clock_t end_tgrowth = std::clock();
//     double elapsed_secs = double(end_tgrowth - begin_tgrowth) / CLOCKS_PER_SEC;
//     BOOST_LOG_SEV(lg, info) << "Growth elapsed secs = " << elapsed_secs;
//     wtime = omp_get_wtime() - wtime;
//     BOOST_LOG_SEV(lg, info) << "Thread number " << omp_get_thread_num()
// 			    << " growth elapsed time " << wtime << " at time = " << ctime;    

//   }//omp master
//   return err;
// }

// ------------------ openmp version ------------------

// WARNING TEST no charging
int NEvo::advance_nocharging(const double ctime) {


  int err;
//   double dtf = 100.0;// time factor multplier
  // with surface growth
  double wsg = static_cast<double>(nm.rs.wsg);
//   // with nucleation
  double wnu = static_cast<double>(nm.rs.wnu);
//   jnucleation[0][chrgs.ineutral] = jnuconst;

//   jnucleation[0][cr.gm.chrgs.maxnegative] = nm.rs.nucleation_rate;

//   // BIMODAL
//   jnucleation[0][cr.gm.chrgs.maxnegative] = nm.rs.nucleation_rate*1E-3;
//   // Warthesen after 10ms no more nucleation, DANGER time splitting and implic
//   if(ctime<1e-2) {
// //  if(tarea<1.8e-2) {
//     wsg = 0.0;
//     wnu = 1.0;//1.0
//     jnucleation[0][cr.gm.chrgs.maxnegative] = nm.rs.nucleation_rate;
//   }
//   // END BIMODAL
  
  jnucleation[0][cr.gm.chrgs.maxnegative] = nm.rs.nucleation_rate;
  // with coagulation
  double wco = static_cast<double>(nm.rs.wco);

  // TODO refactor
  if((nm.rs.wnu == 1) || (nm.rs.wsg == 1)) {
    // Time splitting
    // + nucleation and growth at step dt/2
    err = compute_sgrowth();

    for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
      for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {

        if (pdens[l][q]<0.0){//EPSILON) {
          std::cout << "\n[ee] Negative nanoparticle pdensity (growth + nuc) 1 pdens"
                    << pdens[l][q] << " " << l << " " << q << " .Abort.\n";
          std::cout << "radii = " << cr.gm.vols.radii[l] << "\n charge = " << cr.gm.chrgs.charges[q] << "\n";
          std::terminate();
          pdens[l][q]=0.0;
        }
        ndens[l][q] = (pdens[l][q] + 0.5*nm.tm.ndeltat
                    * (wsg*gsurfacegrowth[l][q]+ wnu*jnucleation[l][q]));

        if (ndens[l][q]<0.0){//EPSILON) {
          std::cout << "\n[ee] Negative nanoparticle ndensity (growth + nuc) 1 ndens";
          std::terminate();
          ndens[l][q]=0.0;
        }
        gsurfacegrowth[l][q] = 0.0;
      }
    }
    // update density
    pdens = ndens;
  }

  clock_t begin_tcoagulation = std::clock();
    
  if(nm.rs.wco == 1) {
    // explicit
    int err = compute_coagulation();

    for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
      for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
        ndens[l][q] = pdens[l][q]
                    + wco * nm.tm.ndeltat
                          * kcoagulation[l][q];
        //DANGER WARNING
        if (ndens[l][q]<0.0/*EPSILONETA*/){//EPSILON) {
  //         std::cout << "\n[ee] Negative nanoparticle ndensity (coagulation)";
          std::cerr << "\n[ee] Negative nanoparticle ndensity (coagulation)";
          //std::terminate();
          ndens[l][q] = 0.0;
        }
        kcoagulation[l][q] = 0.0;
      }
    }
    
    // update density
    pdens = ndens;
    //}
  }


  clock_t end_tcoagulation = std::clock();
  double elapsed_secs = double(end_tcoagulation - begin_tcoagulation) / CLOCKS_PER_SEC;
  BOOST_LOG_SEV(lg, info) << "Coagulation elapsed secs = " << elapsed_secs;

  if((nm.rs.wnu == 1) || (nm.rs.wsg == 1)) {
    // + nucleation and growth at step dt/2
    err = compute_sgrowth();
    for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
      for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {

        if (pdens[l][q]<0.0){//EPSILON) {
          std::cout << "\n[ee] Negative nanoparticle pdensity (growth + nuc) 1 pdens"
                    << pdens[l][q] << " " << l << " " << q << " .Abort.\n";
          std::cout << "radii = " << cr.gm.vols.radii[l] << "\n charge = "
                    << cr.gm.chrgs.charges[q] << "\n";
          std::terminate();
          pdens[l][q]=0.0;
        }
        ndens[l][q] = (pdens[l][q] + 0.5*nm.tm.ndeltat
                    * (wsg*gsurfacegrowth[l][q]+ wnu*jnucleation[l][q]));

        if (ndens[l][q]<0.0){//EPSILON) {
          std::cout << "\n[ee] Negative nanoparticle ndensity (growth + nuc) 1 ndens";
          std::terminate();
          ndens[l][q]=0.0;
        }
        gsurfacegrowth[l][q] = 0.0;
      }
    }
    // update density
    pdens = ndens;
  
  }
  //}
// 
//   // WARNING FIXME
//   total_charge = 0.0;
//   for (unsigned int ll = 0; ll < vols.nsecs; ++ll) {
//     for (unsigned int qq = qind[ll][0]; qq < qind[ll][1]+1; ++qq) {
//       total_charge += chrgs.pivots[qq] * ndens[ll][qq];
//     }
//   }
// 
//   // WARNING time advance
//   ctime += delta_time*dtf;// time factor multiplier
// 
  return 0;
}

// ------------------ openmp version ------------------

// WARNING TEST no charging
int NEvo::advance_nocharging_omp(const double ctime) {


  int err;
//   double dtf = 100.0;// time factor multplier
  // with surface growth
  double wsg = static_cast<double>(nm.rs.wsg);
//   // with nucleation
  double wnu = static_cast<double>(nm.rs.wnu);
//   jnucleation[0][chrgs.ineutral] = jnuconst;

//   jnucleation[0][cr.gm.chrgs.maxnegative] = nm.rs.nucleation_rate;

//   // BIMODAL
//   jnucleation[0][cr.gm.chrgs.maxnegative] = nm.rs.nucleation_rate*1E-3;
//   // Warthesen after 10ms no more nucleation, DANGER time splitting and implic
//   if(ctime<1e-2) {
// //  if(tarea<1.8e-2) {
//     wsg = 0.0;
//     wnu = 1.0;//1.0
//     jnucleation[0][cr.gm.chrgs.maxnegative] = nm.rs.nucleation_rate;
//   }
//   // END BIMODAL
   

  // with coagulation
  double wco = static_cast<double>(nm.rs.wco);

  //#pragma omp single
  //{
  // TODO refactor
  if((nm.rs.wnu == 1) || (nm.rs.wsg == 1)) {
    // Time splitting
    // + nucleation and growth at step dt/2
    if(nm.rs.wsih4==1) {
      jnucleation[0][cr.gm.chrgs.maxnegative] = nm.rs.nucleation_rate;
      err = compute_sgrowth();
      for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
	for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {

	  if (pdens[l][q]<0.0){//EPSILON) {
	    std::cout << "\n[ee] Negative nanoparticle pdensity (growth + nuc) 1 pdens"
		      << pdens[l][q] << " " << l << " " << q << " .Abort.\n";
	    std::cout << "radii = " << cr.gm.vols.radii[l] << "\n charge = " << cr.gm.chrgs.charges[q] << "\n";
	    std::terminate();
	    pdens[l][q]=0.0;
	  }
	  ndens[l][q] = (pdens[l][q] + 0.5*nm.tm.ndeltat
			 * (wsg*gsurfacegrowth[l][q]+wnu*jnucleation[l][q]));

	  // store rates
	  srate2d[l][q] = 0.5*wsg*gsurfacegrowth[l][q];
	  nrate2d[l][q] = 0.5*wnu*jnucleation[l][q];
	
	  if (ndens[l][q]<0.0){//EPSILON) {
	    std::cout << "\n[ee] Negative nanoparticle ndensity (growth + nuc) 1 ndens";
	    std::terminate();
	    ndens[l][q]=0.0;
	  }
	  gsurfacegrowth[l][q] = 0.0;
	}
      }
    }
    else {
      jnucleation[0][cr.gm.chrgs.maxnegative] = nm.rs.nucleation_rate;
      err = compute_sgrowth();
      for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
	for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {

	  if (pdens[l][q]<0.0){//EPSILON) {
	    std::cout << "\n[ee] Negative nanoparticle pdensity (growth + nuc) 1 pdens"
		      << pdens[l][q] << " " << l << " " << q << " .Abort.\n";
	    std::cout << "radii = " << cr.gm.vols.radii[l] << "\n charge = " << cr.gm.chrgs.charges[q] << "\n";
	    std::terminate();
	    pdens[l][q]=0.0;
	  }
	  ndens[l][q] = (pdens[l][q] + 0.5*nm.tm.ndeltat
			 * (wsg*gsurfacegrowth[l][q]+wnu*jnucleation[l][q]));

	  // store rates
	  srate2d[l][q] = 0.5*wsg*gsurfacegrowth[l][q];
	  nrate2d[l][q] = 0.5*wnu*jnucleation[l][q];
	
	  if (ndens[l][q]<0.0){//EPSILON) {
	    std::cout << "\n[ee] Negative nanoparticle ndensity (growth + nuc) 1 ndens";
	    std::terminate();
	    ndens[l][q]=0.0;
	  }
	  gsurfacegrowth[l][q] = 0.0;
	}
      }
    }
    // update density
    pdens = ndens;
  }

  //}
  clock_t begin_tcoagulation = std::clock();
    
  if(nm.rs.wco == 1) {
    // explicit
    int err = compute_coagulation_omp();

    //#pragma omp single
    //{
    for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
      for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
        ndens[l][q] = pdens[l][q]
                    + wco * nm.tm.ndeltat
                          * kcoagulation[l][q];
	crate2d[l][q] = wco * kcoagulation[l][q];
        //DANGER WARNING
        if (ndens[l][q]<0.0/*EPSILONETA*/){//EPSILON) {
          std::cerr << "\n[ee] Negative nanoparticle ndensity (coagulation)";
          //std::terminate();
          ndens[l][q] = 0.0;
        }
        kcoagulation[l][q] = 0.0;
      }
    }
    
    // update density
    pdens = ndens;
    //}
  }


  clock_t end_tcoagulation = std::clock();
  double elapsed_secs = double(end_tcoagulation - begin_tcoagulation) / CLOCKS_PER_SEC;
  BOOST_LOG_SEV(lg, info) << "Coagulation elapsed secs = " << elapsed_secs;

  //#pragma omp single
  //{    
  if((nm.rs.wnu == 1) || (nm.rs.wsg == 1)) {
    // + nucleation and growth at step dt/2
    if(nm.rs.wsih4==1) {
      jnucleation[0][cr.gm.chrgs.maxnegative] = nm.rs.nucleation_rate;
      err = compute_sgrowth();
      for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
	for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {

	  if (pdens[l][q]<0.0){//EPSILON) {
	    std::cout << "\n[ee] Negative nanoparticle pdensity (growth + nuc) 1 pdens"
		      << pdens[l][q] << " " << l << " " << q << " .Abort.\n";
	    std::cout << "radii = " << cr.gm.vols.radii[l] << "\n charge = "
		      << cr.gm.chrgs.charges[q] << "\n";
	    std::terminate();
	    pdens[l][q]=0.0;
	  }
	  ndens[l][q] = (pdens[l][q] + 0.5*nm.tm.ndeltat
			 * (wsg*gsurfacegrowth[l][q]+ wnu*jnucleation[l][q]));

	  // store rates += to add to first split
	  srate2d[l][q] += 0.5*wsg*gsurfacegrowth[l][q];
	  nrate2d[l][q] += 0.5*wnu*jnucleation[l][q];
	
	  if (ndens[l][q]<0.0){//EPSILON) {
	    std::cout << "\n[ee] Negative nanoparticle ndensity (growth + nuc) 1 ndens";
	    std::terminate();
	    ndens[l][q]=0.0;
	  }
	  gsurfacegrowth[l][q] = 0.0;
	}
      }
    }
    else {
      jnucleation[0][cr.gm.chrgs.maxnegative] = nm.rs.nucleation_rate;
      err = compute_sgrowth();
      for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
	for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {

	  if (pdens[l][q]<0.0){//EPSILON) {
	    std::cout << "\n[ee] Negative nanoparticle pdensity (growth + nuc) 1 pdens"
		      << pdens[l][q] << " " << l << " " << q << " .Abort.\n";
	    std::cout << "radii = " << cr.gm.vols.radii[l] << "\n charge = "
		      << cr.gm.chrgs.charges[q] << "\n";
	    std::terminate();
	    pdens[l][q]=0.0;
	  }
	  ndens[l][q] = (pdens[l][q] + 0.5*nm.tm.ndeltat
			 * (wsg*gsurfacegrowth[l][q]+ wnu*jnucleation[l][q]));

	  // store rates += to add to first split
	  srate2d[l][q] += 0.5*wsg*gsurfacegrowth[l][q];
	  nrate2d[l][q] += 0.5*wnu*jnucleation[l][q];
	
	  if (ndens[l][q]<0.0){//EPSILON) {
	    std::cout << "\n[ee] Negative nanoparticle ndensity (growth + nuc) 1 ndens";
	    std::terminate();
	    ndens[l][q]=0.0;
	  }
	  gsurfacegrowth[l][q] = 0.0;
	}
      }
    }
    // update density
    pdens = ndens;
  
  }
  //}
// 
//   // WARNING FIXME
//   total_charge = 0.0;
//   for (unsigned int ll = 0; ll < vols.nsecs; ++ll) {
//     for (unsigned int qq = qind[ll][0]; qq < qind[ll][1]+1; ++qq) {
//       total_charge += chrgs.pivots[qq] * ndens[ll][qq];
//     }
//   }
// 
//   // WARNING time advance
//   ctime += delta_time*dtf;// time factor multiplier
// 
  return 0;
}


// ------------------ openmp version ------------------

// WARNING TEST no charging
int NEvo::advance_nocharging_ompadp(const double ctime) {

  int err;
//   double dtf = 100.0;// time factor multplier
  // with surface growth
  double wsg = static_cast<double>(nm.rs.wsg);
//   // with nucleation
  double wnu = static_cast<double>(nm.rs.wnu);
//   jnucleation[0][chrgs.ineutral] = jnuconst;

//   jnucleation[0][cr.gm.chrgs.maxnegative] = nm.rs.nucleation_rate;

//   // BIMODAL
//   jnucleation[0][cr.gm.chrgs.maxnegative] = nm.rs.nucleation_rate*1E-3;
//   // Warthesen after 10ms no more nucleation, DANGER time splitting and implic
//   if(ctime<1e-2) {
// //  if(tarea<1.8e-2) {
//     wsg = 0.0;
//     wnu = 1.0;//1.0
//     jnucleation[0][cr.gm.chrgs.maxnegative] = nm.rs.nucleation_rate;
//   }
//   // END BIMODAL

  jnucleation[0][cr.gm.chrgs.maxnegative] = nm.rs.nucleation_rate;

  // with coagulation
  double wco = static_cast<double>(nm.rs.wco);

  //#pragma omp single
  //{
  // TODO refactor
  // if((nm.rs.wnu == 1) || (nm.rs.wsg == 1)) {
  if((wnu == 1) || (wsg == 1)) {
    compute_split_sgnucleation();
  }

  //}
  double begin_tcoagulation = omp_get_wtime();
    
  if(nm.rs.wco == 1) {
    // explicit
    int err = compute_coagulation_omp();

    //#pragma omp single
    //{
    //#pragma omp parallel
    //{
      //#pragma omp for nowait collapse(2)
    //#pragma omp parallel for
    for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
      //#pragma omp parallel for
      for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
        ndens[l][q] = pdens[l][q]
                    + nm.rs.wco * nm.tm.ndeltat
                          * kcoagulation[l][q];
	crate2d[l][q] = nm.rs.wco * kcoagulation[l][q];
        //DANGER WARNING
        if (ndens[l][q]<0.0/*EPSILONETA*/){//EPSILON) {
          //std::cerr << "\n[ee] Negative nanoparticle ndensity (coagulation)";
	  //BOOST_LOG_SEV(lg, info) << "\n[ee] Negative nanoparticle ndensity (coagulation) : "
	  //			  << ndens[l][q] << " , " << pdens[l][q];
	  //std::terminate();
	  //std::cerr << "\n[ee] Diameter " << cr.gm.vols.diameters[l]*2e9;
	  //std::cerr << "\n[ee] Charge " << cr.gm.chrgs.charges[q];
	  //std::cerr << "\n[ee] pdens " << pdens[l][q];
	  //std::cerr << "\n[ee] ndens " << ndens[l][q];
	  //std::cerr << "\n[ee]";
	  // if (ndens[l][q]<-1.0e-2) {
	  //   nm.tm.ndeltat *= 0.75;
	  //   BOOST_LOG_SEV(lg, info) << "\n[ee] New dt = ";
	  // }
	  if (ndens[l][q]<-1.0e-1) {
	    BOOST_LOG_SEV(lg, info) << "Negative nanoparticle ndensity (coagulation)";
	    //nm.tm.ndeltat *= 0.9;
	    BOOST_LOG_SEV(lg, info) << "New dt = " << nm.tm.ndeltat;
	    if (nm.tm.ndeltat < nm.tm.qdeltat) {
	      //std::terminate();
	    }
	  }
	  ndens[l][q] = 0.0;
	  if (pdens[l][q] > 0.0) {
	    ndens[l][q] = pdens[l][q];
	  }
          
	  
        }
        kcoagulation[l][q] = 0.0;
      }// for q 
    }// for l
    //    }// omp parallel
    // update density
    //#pragma omp barrier
    pdens = ndens;
    //}
  }

  double end_tcoagulation = omp_get_wtime();
  double elapsed_secs = end_tcoagulation - begin_tcoagulation;
  BOOST_LOG_SEV(lg, info) << "Coagulation elapsed secs = " << elapsed_secs;

  //#pragma omp single
  //{    
  if((nm.rs.wnu == 1) || (nm.rs.wsg == 1)) {
    compute_split_sgnucleation();  
  }
  //}
  return 0;
}

// ------------------ adaptive ------------------

// WARNING TEST no charging
int NEvo::advance_nocharging_adapt(const double ctime) {

  int err;
//   double dtf = 100.0;// time factor multplier
  // with surface growth
  double wsg = static_cast<double>(nm.rs.wsg);
//   // with nucleation
  double wnu = static_cast<double>(nm.rs.wnu);
//   jnucleation[0][chrgs.ineutral] = jnuconst;

//   jnucleation[0][cr.gm.chrgs.maxnegative] = nm.rs.nucleation_rate;

//   // BIMODAL
//   jnucleation[0][cr.gm.chrgs.maxnegative] = nm.rs.nucleation_rate*1E-3;
//   // Warthesen after 10ms no more nucleation, DANGER time splitting and implic
//   if(ctime<1e-2) {
// //  if(tarea<1.8e-2) {
//     wsg = 0.0;
//     wnu = 1.0;//1.0
//     jnucleation[0][cr.gm.chrgs.maxnegative] = nm.rs.nucleation_rate;
//   }
//   // END BIMODAL

  jnucleation[0][cr.gm.chrgs.maxnegative] = nm.rs.nucleation_rate;

  // with coagulation
  double wco = static_cast<double>(nm.rs.wco);

  //#pragma omp single
  //{
  // TODO refactor
  // if((nm.rs.wnu == 1) || (nm.rs.wsg == 1)) {
  if((wnu == 1) || (wsg == 1)) {
    compute_split_sgnucleation();
  }

  //}
  double begin_tcoagulation = omp_get_wtime();
    
  if(nm.rs.wco == 1) {
    // explicit
    int err = compute_coagulation_omp();

    //#pragma omp single
    //{
    //#pragma omp parallel
    //{
      //#pragma omp for nowait collapse(2)
    //#pragma omp parallel for
    for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
      //#pragma omp parallel for
      for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
        ndens[l][q] = pdens[l][q]
                    + nm.rs.wco * nm.tm.ndeltat
                          * kcoagulation[l][q];
	crate2d[l][q] = nm.rs.wco * kcoagulation[l][q];
        //DANGER WARNING
        if (ndens[l][q]<0.0/*EPSILONETA*/){//EPSILON) {
          //std::cerr << "\n[ee] Negative nanoparticle ndensity (coagulation)";
	  //BOOST_LOG_SEV(lg, info) << "\n[ee] Negative nanoparticle ndensity (coagulation) : "
	  //			  << ndens[l][q] << " , " << pdens[l][q];
	  //std::terminate();
	  //std::cerr << "\n[ee] Diameter " << cr.gm.vols.diameters[l]*2e9;
	  //std::cerr << "\n[ee] Charge " << cr.gm.chrgs.charges[q];
	  //std::cerr << "\n[ee] pdens " << pdens[l][q];
	  //std::cerr << "\n[ee] ndens " << ndens[l][q];
	  //std::cerr << "\n[ee]";
	  // if (ndens[l][q]<-1.0e-2) {
	  //   nm.tm.ndeltat *= 0.75;
	  //   BOOST_LOG_SEV(lg, info) << "\n[ee] New dt = ";
	  // }
	  if (ndens[l][q]<-1.0e-1) {
	    BOOST_LOG_SEV(lg, info) << "Negative nanoparticle ndensity (coagulation)";
	    //nm.tm.ndeltat *= 0.9;
	    BOOST_LOG_SEV(lg, info) << "New dt = " << nm.tm.ndeltat;
	    if (nm.tm.ndeltat < nm.tm.qdeltat) {
	      //std::terminate();
	    }
	  }
	  ndens[l][q] = 0.0;
	  if (pdens[l][q] > 0.0) {
	    ndens[l][q] = pdens[l][q];
	  }
          
	  
        }
        kcoagulation[l][q] = 0.0;
      }// for q 
    }// for l
    //    }// omp parallel
    // update density
    //#pragma omp barrier
    pdens = ndens;
    //}
  }

  double end_tcoagulation = omp_get_wtime();
  double elapsed_secs = end_tcoagulation - begin_tcoagulation;
  BOOST_LOG_SEV(lg, info) << "Coagulation elapsed secs = " << elapsed_secs;

  //#pragma omp single
  //{    
  if((nm.rs.wnu == 1) || (nm.rs.wsg == 1)) {
    compute_split_sgnucleation();  
  }
  //}
  return 0;
}



// ========================= SiH4 =========================
int NEvo::advance_nocharging_ompadpsih4(const double ctime) {

  int err;

  darray new_mom(0.0, 3);
  darray old_mom(0.0, 3);
  
  // with surface growth
  double wsg = static_cast<double>(nm.rs.wsg);
  // with nucleation
  double wnu = static_cast<double>(nm.rs.wnu);

  jnucleation[0][cr.gm.chrgs.maxnegative] = nm.rs.nucleation_rate;

  // with coagulation
  double wco = static_cast<double>(nm.rs.wco);

  // compute_moments(pdens, old_mom);
  
  if((wnu == 1) || (wsg == 1)) {
    compute_split_sgnucleationsih4();
  }

  // compute_moments(pdens, new_mom);
  
  // BOOST_LOG_SEV(lg, info) << "number rel error " << ctime << '\t' << relative_error<double>(old_mom[0], new_mom[0]);
  // BOOST_LOG_SEV(lg, info) << "number abs error " << ctime << '\t' << absolute_error<double>(old_mom[0], new_mom[0]);
  // BOOST_LOG_SEV(lg, info) << "volume rel error " << ctime << '\t' << relative_error<double>(old_mom[1], new_mom[1]);
  // BOOST_LOG_SEV(lg, info) << "volume abs error " << ctime << '\t' << absolute_error<double>(old_mom[1], new_mom[1]);
  // BOOST_LOG_SEV(lg, info) << "charge rel error " << ctime << '\t' << relative_error<double>(old_mom[2], new_mom[2]);
  // BOOST_LOG_SEV(lg, info) << "charge abs error " << ctime << '\t' << absolute_error<double>(old_mom[2], new_mom[2]);
    
  //}
  double begin_tcoagulation = omp_get_wtime();
    
  if(nm.rs.wco == 1) {

    compute_moments(pdens, old_mom);

    // explicit
    int err = compute_coagulation_omp();

    //#pragma omp single
    //{
    //#pragma omp parallel
    //{
      //#pragma omp for nowait collapse(2)
    //#pragma omp parallel for
    for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
      //#pragma omp parallel for
      for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
        ndens[l][q] = pdens[l][q]
                    + nm.rs.wco * nm.tm.ndeltat
                          * kcoagulation[l][q];
	crate2d[l][q] = nm.rs.wco * kcoagulation[l][q];
        //DANGER WARNING
        if (ndens[l][q]<0.0/*EPSILONETA*/){//EPSILON) {
          //std::cerr << "\n[ee] Negative nanoparticle ndensity (coagulation)";
	  //BOOST_LOG_SEV(lg, info) << "\n[ee] Negative nanoparticle ndensity (coagulation) : "
	  //			  << ndens[l][q] << " , " << pdens[l][q];
	  //std::terminate();
	  //std::cerr << "\n[ee] Diameter " << cr.gm.vols.diameters[l]*2e9;
	  //std::cerr << "\n[ee] Charge " << cr.gm.chrgs.charges[q];
	  //std::cerr << "\n[ee] pdens " << pdens[l][q];
	  //std::cerr << "\n[ee] ndens " << ndens[l][q];
	  //std::cerr << "\n[ee]";
	  // if (ndens[l][q]<-1.0e-2) {
	  //   nm.tm.ndeltat *= 0.75;
	  //   BOOST_LOG_SEV(lg, info) << "\n[ee] New dt = ";
	  // }
	  if (ndens[l][q]<-1.0e-1) {
	    BOOST_LOG_SEV(lg, info) << "Negative nanoparticle ndensity (coagulation)";
	    nm.tm.ndeltat *= 0.9;
	    BOOST_LOG_SEV(lg, info) << "New dt = " << nm.tm.ndeltat;
	    if (nm.tm.ndeltat < nm.tm.qdeltat) {
	      std::terminate();
	    }
	  }
	  ndens[l][q] = 0.0;
	  if (pdens[l][q] > 0.0) {
	    ndens[l][q] = pdens[l][q];
	  }
          
	  
        }
        kcoagulation[l][q] = 0.0;
      }// for q 
    }// for l
    //    }// omp parallel
    // update density
    //#pragma omp barrier

    // compute_moments(ndens, new_mom);

    // BOOST_LOG_SEV(lg, info) << "number rel error " << ctime << '\t' << relative_error<double>(old_mom[0], new_mom[0]);
    // BOOST_LOG_SEV(lg, info) << "number abs error " << ctime << '\t' << absolute_error<double>(old_mom[0], new_mom[0]);
    // BOOST_LOG_SEV(lg, info) << "volume rel error " << ctime << '\t' << relative_error<double>(old_mom[1], new_mom[1]);
    // BOOST_LOG_SEV(lg, info) << "volume abs error " << ctime << '\t' << absolute_error<double>(old_mom[1], new_mom[1]);
    // BOOST_LOG_SEV(lg, info) << "charge rel error " << ctime << '\t' << relative_error<double>(old_mom[2], new_mom[2]);
    // BOOST_LOG_SEV(lg, info) << "charge abs error " << ctime << '\t' << absolute_error<double>(old_mom[2], new_mom[2]);

    // if(is_close(new_mom[1], old_mom[1], 1e-2)==false) {
    //   std::cerr << "\nError in coagulation, moments differ\n";
    //   BOOST_LOG_SEV(lg, info) << "Error in coagulation, moments differ";
    //   BOOST_LOG_SEV(lg, info) << "pdens N = " << old_mom[0];
    //   BOOST_LOG_SEV(lg, info) << "ndens N = " << new_mom[0];
    //   BOOST_LOG_SEV(lg, info) << "cn0 " << abs(new_mom[0]-old_mom[0]);
    //   BOOST_LOG_SEV(lg, info) << "pdens V = " << old_mom[1];
    //   BOOST_LOG_SEV(lg, info) << "ndens V = " << new_mom[1];
    //   BOOST_LOG_SEV(lg, info) << "mv0 " << abs(new_mom[1]-old_mom[1]);
    //   std::terminate();
    // }
    
    pdens = ndens;
    //}
  }

  double end_tcoagulation = omp_get_wtime();
  double elapsed_secs = end_tcoagulation - begin_tcoagulation;
  BOOST_LOG_SEV(lg, info) << "Coagulation elapsed secs = " << elapsed_secs;

  //#pragma omp single
  //{    
  if((nm.rs.wnu == 1) || (nm.rs.wsg == 1)) {
    compute_split_sgnucleationsih4();  
  }
  //}
  return 0;
}


// ========================= SiH4 =========================
int NEvo::advance_nosplit_ompadpsih4(const double ctime) {

  int err;

  darray new_mom(0.0, 3);
  darray old_mom(0.0, 3);
  
  // with surface growth
  double wsg = static_cast<double>(nm.rs.wsg);
  // with nucleation
  double wnu = static_cast<double>(nm.rs.wnu);

  //jnucleation[0][cr.gm.chrgs.maxnegative] = nm.rs.nucleation_rate;

  // with coagulation
  double wco = static_cast<double>(nm.rs.wco);

  // compute_moments(pdens, old_mom);
  
  // if((wnu == 1) || (wsg == 1)) {
  //   compute_split_sgnucleationsih4();
  // }

  // compute_moments(pdens, new_mom);
  
  // BOOST_LOG_SEV(lg, info) << "number rel error " << ctime << '\t' << relative_error<double>(old_mom[0], new_mom[0]);
  // BOOST_LOG_SEV(lg, info) << "number abs error " << ctime << '\t' << absolute_error<double>(old_mom[0], new_mom[0]);
  // BOOST_LOG_SEV(lg, info) << "volume rel error " << ctime << '\t' << relative_error<double>(old_mom[1], new_mom[1]);
  // BOOST_LOG_SEV(lg, info) << "volume abs error " << ctime << '\t' << absolute_error<double>(old_mom[1], new_mom[1]);
  // BOOST_LOG_SEV(lg, info) << "charge rel error " << ctime << '\t' << relative_error<double>(old_mom[2], new_mom[2]);
  // BOOST_LOG_SEV(lg, info) << "charge abs error " << ctime << '\t' << absolute_error<double>(old_mom[2], new_mom[2]);
    
  //}
  double begin_tcoagulation = omp_get_wtime();
    
  if(nm.rs.wco == 1) {

    compute_moments(pdens, old_mom);

    // explicit
    int err = compute_coagulation_omp();

    err = compute_nosplit_sgnucleationsih4();

    //#pragma omp single
    //{
    //#pragma omp parallel
    //{
      //#pragma omp for nowait collapse(2)
    //#pragma omp parallel for
    darray radii = cr.gm.vols.radii;
    darray charges = cr.gm.chrgs.charges;


    // // compute total rate x dt
    // for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
    //   //#pragma omp parallel for
    //   for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
    //     trate2d[l][q] = (wco*kcoagulation[l][q]+wsg*gsurfacegrowth[l][q]+wnu*jnucleation[l][q]);
    // 	double fract = std::abs(pdens[l][q]/(trate2d[l][q]));
    // 	if ((nm.tm.ndeltat*std::abs(trate2d[l][q]) > pdens[l][q]) && (pdens[l][q] > 1)
    // 	    && (trate2d[l][q] < 0)) {
    // 	  dt2d[l][q] = std::pow(10, std::floor(std::log10(fract)));
    // 	  BOOST_LOG_SEV(lg, info) << "DT l,q, rl, qq\t" << l << '\t' << q << '\t'
    // 				  << radii[l]*1e9 << '\t' << charges[q];
    // 	  BOOST_LOG_SEV(lg, info) << "fract\t" << fract;
    // 	  BOOST_LOG_SEV(lg, info) << "DT\t" << dt2d[l][q];
    // 	  BOOST_LOG_SEV(lg, info) << "trate\t" << trate2d[l][q];
    // 	  BOOST_LOG_SEV(lg, info) << "pdens\t" << pdens[l][q];
    // 	  BOOST_LOG_SEV(lg, info) << "Trate\t" << nm.tm.ndeltat*trate2d[l][q];
    // 	  BOOST_LOG_SEV(lg, info) << "newTrate\t" << dt2d[l][q]*trate2d[l][q];
    // 	  BOOST_LOG_SEV(lg, info) << "newdens\t" << pdens[l][q]+ dt2d[l][q]*trate2d[l][q];
    // 	  BOOST_LOG_SEV(lg, info) << "steps\t" <<  nm.tm.ndeltat/dt2d[l][q];
    // 	}
    // 	else {
    // 	  dt2d[l][q] = nm.tm.ndeltat;
    // 	}
    //   }
    // }

    // BOOST_LOG_SEV(lg, info) << "ctiming " <<  ctime;

    // //INIT:
    // // compute total rate x dt
    // for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
    //   //#pragma omp parallel for
    //   for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
    // 	size_t steps = static_cast<size_t>(nm.tm.ndeltat/dt2d[l][q]);
    // 	for(size_t tr = 0; tr < (steps<1?1:steps); ++tr) {
    // 	  ndens[l][q] = pdens[l][q] + dt2d[l][q] * trate2d[l][q];
	  
    // 	  // crate2d[l][q] = wco * kcoagulation[l][q];
    // 	  // //DANGER WARNING
    // 	  // if (ndens[l][q]<0.0/*EPSILONETA*/){//EPSILON) {
    // 	  //   if (ndens[l][q]<-1.0e-1) {
    // 	  //     BOOST_LOG_SEV(lg, info) << "Negative nanoparticle ndensity (coagulation)";
    // 	  //     BOOST_LOG_SEV(lg, info) << ndens[l][q];
    // 	  //     BOOST_LOG_SEV(lg, info) << pdens[l][q];
    // 	  //     BOOST_LOG_SEV(lg, info) << wco*kcoagulation[l][q];
    // 	  //     BOOST_LOG_SEV(lg, info) << wsg*gsurfacegrowth[l][q];
    // 	  //     BOOST_LOG_SEV(lg, info) << wnu*jnucleation[l][q];
    // 	  //     BOOST_LOG_SEV(lg, info) << wco*kcoagulation[l][q]+wsg*gsurfacegrowth[l][q]+wnu*jnucleation[l][q];
    // 	  //     ndens[l][q] = 0.0;
    // 	  //     //std::terminate();
    // 	  //     //nm.tm.ndeltat *= 0.9;
    // 	  //     BOOST_LOG_SEV(lg, info) << "t, New dt = " << ctime << '\t' << nm.tm.ndeltat;
    // 	  //     BOOST_LOG_SEV(lg, info) << "l,q, rl, qq\t" << l << '\t' << q << '\t'
    // 	  // 			      << radii[l]*1e9 << '\t' << charges[q];
    // 	  //     if (nm.tm.ndeltat < nm.tm.qdeltat) {
    // 	  // 	std::terminate();
    // 	  //     }
    // 	  //     std::terminate();
    // 	  //     //goto INIT;
    // 	  //   }
    // 	  //   ndens[l][q] = 0.0;
    // 	  //   if (pdens[l][q] > 0.0) {
    // 	  //     ndens[l][q] = pdens[l][q];
    // 	  //   }
    // 	  // }
	  
    // 	  pdens[l][q] = ndens[l][q];
    // 	}
    // 	kcoagulation[l][q] = 0.0;
    //   }// for q 
    // }// for l


    for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
      //#pragma omp parallel for
      for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
        ndens[l][q] = pdens[l][q]
  	  + nm.tm.ndeltat
  	  * (wco*kcoagulation[l][q]+wsg*gsurfacegrowth[l][q]+wnu*jnucleation[l][q]);

	// // semi implicit
	// ndens[l][q] = (pdens[l][q]
  	//   + nm.tm.ndeltat
	// 	       * (wco*birth_vector[l][q]+wsg*gsurfacegrowth[l][q]+wnu*jnucleation[l][q]))
	//   /(1.0+nm.tm.ndeltat*death_vector[l][q]);

  	crate2d[l][q] = wco * kcoagulation[l][q];
        //DANGER WARNING
        if (ndens[l][q]<0.0/*EPSILONETA*/){//EPSILON) {
  	  if (ndens[l][q]<-1.0e-1) {
	    BOOST_LOG_SEV(lg, info) << "\n[ee] Negative nanoparticle ndensity (coagulation) : "
				    << ndens[l][q] << " , " << pdens[l][q];

	    BOOST_LOG_SEV(lg, info) << "\n[ee] Diameter " << cr.gm.vols.diameters[l]*2e9;
	    BOOST_LOG_SEV(lg, info) << "\n[ee] Charge " << cr.gm.chrgs.charges[q];
   	    BOOST_LOG_SEV(lg, info) << ndens[l][q];
  	    BOOST_LOG_SEV(lg, info) << pdens[l][q];
    	    BOOST_LOG_SEV(lg, info) << wco*kcoagulation[l][q];
  	    BOOST_LOG_SEV(lg, info) << wsg*gsurfacegrowth[l][q];
  	    BOOST_LOG_SEV(lg, info) << wnu*jnucleation[l][q];
  	    BOOST_LOG_SEV(lg, info) << wco*kcoagulation[l][q]+wsg*gsurfacegrowth[l][q]+wnu*jnucleation[l][q];
  	    if (nm.tm.ndeltat < nm.tm.qdeltat) {
  	      std::terminate();
  	    }
	    return -1;
  	  }
  	  ndens[l][q] = 0.0;
  	  if (pdens[l][q] > 0.0) {
  	    ndens[l][q] = pdens[l][q];
  	  }
        }
        kcoagulation[l][q] = 0.0;
      }// for q 
    }// for l
    
  // INIT:
  //   for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
  //     //#pragma omp parallel for
  //     for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
  //       ndens[l][q] = pdens[l][q]// / 1e10
  // 	  + nm.tm.ndeltat// / 1e10
  // 	  * (wco*kcoagulation[l][q]+wsg*gsurfacegrowth[l][q]+wnu*jnucleation[l][q]);
  // 	// for (unsigned int ll = 0; ll < cr.gm.vols.nsections; ++ll) {
  // 	//   for (unsigned int qq = 0; qq < cr.gm.chrgs.nsections; ++qq) {
  // 	//     ndens[ll][qq] *= 1e10;
  // 	//   }
  // 	// }
  // 	crate2d[l][q] = wco * kcoagulation[l][q];
  //       //DANGER WARNING
  //       if (ndens[l][q]<0.0/*EPSILONETA*/){//EPSILON) {
  //         //std::cerr << "\n[ee] Negative nanoparticle ndensity (coagulation)";
  // 	  //BOOST_LOG_SEV(lg, info) << "\n[ee] Negative nanoparticle ndensity (coagulation) : "
  // 	  //			  << ndens[l][q] << " , " << pdens[l][q];
  // 	  //std::terminate();
  // 	  //std::cerr << "\n[ee] Diameter " << cr.gm.vols.diameters[l]*2e9;
  // 	  //std::cerr << "\n[ee] Charge " << cr.gm.chrgs.charges[q];
  // 	  //std::cerr << "\n[ee] pdens " << pdens[l][q];
  // 	  //std::cerr << "\n[ee] ndens " << ndens[l][q];
  // 	  //std::cerr << "\n[ee]";
  // 	  // if (ndens[l][q]<-1.0e-2) {
  // 	  //   nm.tm.ndeltat *= 0.75;
  // 	  //   BOOST_LOG_SEV(lg, info) << "\n[ee] New dt = ";
  // 	  // }
  // 	  if (ndens[l][q]<-1.0e-1) {
  // 	    BOOST_LOG_SEV(lg, info) << "Negative nanoparticle ndensity (coagulation)";
  // 	    BOOST_LOG_SEV(lg, info) << ndens[l][q];
  // 	    BOOST_LOG_SEV(lg, info) << pdens[l][q];
  //   	    BOOST_LOG_SEV(lg, info) << wco*kcoagulation[l][q];
  // 	    BOOST_LOG_SEV(lg, info) << wsg*gsurfacegrowth[l][q];
  // 	    BOOST_LOG_SEV(lg, info) << wnu*jnucleation[l][q];
  // 	    BOOST_LOG_SEV(lg, info) << wco*kcoagulation[l][q]+wsg*gsurfacegrowth[l][q]+wnu*jnucleation[l][q];
  // 	    ndens[l][q] = 0.0;
  // 	    //std::terminate();
  // 	    nm.tm.ndeltat *= 0.95;
  // 	    BOOST_LOG_SEV(lg, info) << "t, New dt = " << ctime << '\t' << nm.tm.ndeltat;
  // 	    BOOST_LOG_SEV(lg, info) << "l,q, rl, qq\t" << l << '\t' << q << '\t'
  // 				    << radii[l]*1e9 << '\t' << charges[q];
  // 	    if (nm.tm.ndeltat < nm.tm.qdeltat) {
  // 	      std::terminate();
  // 	    }
  // 	    goto INIT;
  // 	  }
  // 	  ndens[l][q] = 0.0;
  // 	  if (pdens[l][q] > 0.0) {
  // 	    ndens[l][q] = pdens[l][q];
  // 	  }
  //       }
  //       kcoagulation[l][q] = 0.0;
  //     }// for q 
  //   }// for l
       // }// omp parallel
    // update density
    //#pragma omp barrier

    // compute_moments(ndens, new_mom);

    // // BOOST_LOG_SEV(lg, info) << "number rel error " << ctime << '\t' << relative_error<double>(old_mom[0], new_mom[0]);
    // // BOOST_LOG_SEV(lg, info) << "number abs error " << ctime << '\t' << absolute_error<double>(old_mom[0], new_mom[0]);
    // // BOOST_LOG_SEV(lg, info) << "volume rel error " << ctime << '\t' << relative_error<double>(old_mom[1], new_mom[1]);
    // // BOOST_LOG_SEV(lg, info) << "volume abs error " << ctime << '\t' << absolute_error<double>(old_mom[1], new_mom[1]);
    // // BOOST_LOG_SEV(lg, info) << "charge rel error " << ctime << '\t' << relative_error<double>(old_mom[2], new_mom[2]);
    // // BOOST_LOG_SEV(lg, info) << "charge abs error " << ctime << '\t' << absolute_error<double>(old_mom[2], new_mom[2]);

    // if(is_close(new_mom[1], old_mom[1], 1e-2)==false) {
    //   std::cerr << "\nError in coagulation, moments differ\n";
    //   BOOST_LOG_SEV(lg, info) << "Error in coagulation, moments differ";
    //   BOOST_LOG_SEV(lg, info) << "pdens N = " << old_mom[0];
    //   BOOST_LOG_SEV(lg, info) << "ndens N = " << new_mom[0];
    //   BOOST_LOG_SEV(lg, info) << "cn0 " << abs(new_mom[0]-old_mom[0]);
    //   BOOST_LOG_SEV(lg, info) << "pdens V = " << old_mom[1];
    //   BOOST_LOG_SEV(lg, info) << "ndens V = " << new_mom[1];
    //   BOOST_LOG_SEV(lg, info) << "mv0 " << abs(new_mom[1]-old_mom[1]);
    //   //std::terminate();
    // }
    
    pdens = ndens;
    //}
  }

  // double end_tcoagulation = omp_get_wtime();
  // double elapsed_secs = end_tcoagulation - begin_tcoagulation;
  // BOOST_LOG_SEV(lg, info) << "Coagulation elapsed secs = " << elapsed_secs
  // 			  << "\t t=" << ctime;

  // //#pragma omp single
  // //{    
  // if((nm.rs.wnu == 1) || (nm.rs.wsg == 1)) {
  //   compute_split_sgnucleationsih4();  
  // }
  // //}
  return 0;
}


// full adaptive no split growth

// ========================= SiH4 =========================
int NEvo::advance_nosplit_sih4_adapt(const double ctime) {

  int err;

  darray new_mom(0.0, 3);
  darray old_mom(0.0, 3);
  
  // with surface growth
  double wsg = static_cast<double>(nm.rs.wsg);
  // with nucleation
  double wnu = static_cast<double>(nm.rs.wnu);

  //jnucleation[0][cr.gm.chrgs.maxnegative] = nm.rs.nucleation_rate;

  // with coagulation
  double wco = static_cast<double>(nm.rs.wco);

  double begin_tcoagulation = omp_get_wtime();
    
  if(nm.rs.wco == 1) {

    compute_moments(pdens, old_mom);

    // explicit
    int err = compute_coagulation_omp();

    err = compute_nosplit_sgnucleationsih4();

    darray radii = cr.gm.vols.radii;
    darray charges = cr.gm.chrgs.charges;


    for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
      //#pragma omp parallel for
      for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
        ndens[l][q] = pdens[l][q]// / 1e10
  	  + nm.tm.ndeltat// / 1e10
  	  * (wco*kcoagulation[l][q]+wsg*gsurfacegrowth[l][q]+wnu*jnucleation[l][q]);

  	crate2d[l][q] = wco * kcoagulation[l][q];
        //DANGER WARNING
        if (ndens[l][q]<0.0/*EPSILONETA*/){//EPSILON) {
  	  if (ndens[l][q]<-1.0e-1) {
	    BOOST_LOG_SEV(lg, info) << "\n[ee] Negative nanoparticle ndensity (coagulation) : "
				    << ndens[l][q] << " , " << pdens[l][q];

	    BOOST_LOG_SEV(lg, info) << "\n[ee] Diameter " << cr.gm.vols.diameters[l]*2e9;
	    BOOST_LOG_SEV(lg, info) << "\n[ee] Charge " << cr.gm.chrgs.charges[q];
   	    BOOST_LOG_SEV(lg, info) << ndens[l][q];
  	    BOOST_LOG_SEV(lg, info) << pdens[l][q];
    	    BOOST_LOG_SEV(lg, info) << wco*kcoagulation[l][q];
  	    BOOST_LOG_SEV(lg, info) << wsg*gsurfacegrowth[l][q];
  	    BOOST_LOG_SEV(lg, info) << wnu*jnucleation[l][q];
  	    BOOST_LOG_SEV(lg, info) << wco*kcoagulation[l][q]+wsg*gsurfacegrowth[l][q]+wnu*jnucleation[l][q];
  	    if (nm.tm.ndeltat < nm.tm.qdeltat) {
  	      std::terminate();
  	    }
	    return -1;
  	  }
  	  ndens[l][q] = 0.0;
  	  if (pdens[l][q] > 0.0) {
  	    ndens[l][q] = pdens[l][q];
  	  }
        }
        kcoagulation[l][q] = 0.0;
      }// for q 
    }// for l    
    
    pdens = ndens;

  }

  // double end_tcoagulation = omp_get_wtime();
  // double elapsed_secs = end_tcoagulation - begin_tcoagulation;
  // BOOST_LOG_SEV(lg, info) << "Coagulation elapsed secs = " << elapsed_secs
  // 			  << "\t t=" << ctime;

  // //#pragma omp single
  // //{    
  // if((nm.rs.wnu == 1) || (nm.rs.wsg == 1)) {
  //   compute_split_sgnucleationsih4();  
  // }
  // //}
  return 0;
}

// ------------------ adaptive sgrowth ----------------

int NEvo::compute_split_sgnucleation() {
  // Time splitting
  // + nucleation and growth at step dt/2
    
  double wsg = static_cast<double>(nm.rs.wsg);
  double wnu = static_cast<double>(nm.rs.wnu);

  unsigned int max_iter = 11;
  unsigned int inner_iter = 0;
  
  // delta time for particle growth (split)
  double dt_growth = 0.5 * nm.tm.ndeltat;
  
  // delta time for sg and nucleation
  double dt_sgn = dt_growth;
  double dt_inner = dt_growth;

  bool negative_density = false;

  int err;
 
  do {
    inner_iter += 1;
    negative_density = false;
    pdens_aux = pdens;
    ndens_aux = zero2d;

    srate2d_aux = zero2d;
    nrate2d_aux = zero2d;
    dt_inner = 0.0;
    do {
      dt_inner += dt_sgn;
      //std::cout << "\n[ee] dt_sgn " << dt_sgn << " inner: " << dt_inner << " iter: " << inner_iter;
      err = compute_sgrowth_adp();
      
      for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
	for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
	  ndens_aux[l][q] = (pdens_aux[l][q] + dt_sgn
			  * (wsg*gsurfacegrowth[l][q]+wnu*jnucleation[l][q]));
	
	  // // store rates
	  // if (first_split) {
	  //   srate2d[l][q] = 0.5*wsg*gsurfacegrowth[l][q];
	  //   nrate2d[l][q] = 0.5*wnu*jnucleation[l][q];
	  // }
	  // else {
	  //   srate2d[l][q] += 0.5*wsg*gsurfacegrowth[l][q];
	  //   nrate2d[l][q] += 0.5*wnu*jnucleation[l][q];
	  // }
	  srate2d_aux[l][q] += 0.5*wsg*gsurfacegrowth[l][q];
	  //nrate2d_aux[l][q] += wnu*jnucleation[l][q];

	  if (ndens_aux[l][q]<0.0 && fabs(ndens_aux[l][q])<1.0) {
	    ndens_aux[l][q] = 0.0;
	  }
	  if (ndens_aux[l][q]<0.0){//EPSILON) {
	    std::cout << "\n[ee] Negative nanoparticle ndensity (growth + nuc) 1 ndens ("
		      << l << ", " << q << ") dt_sgn = " << dt_sgn
		      << " ndens[l][q] = " << ndens_aux[l][q];
	    //std::terminate();
	    negative_density = true;
	    dt_inner = 0.0;
	    break;
	  }
	  //gsurfacegrowth[l][q] = 0.0;
	}
	if (negative_density) {
	  std::cout << "\n[ee] Break 1 : " << inner_iter;
	  break;
	}
      }
      if (negative_density) {
	std::cout << "\n[ee] Break 2 : "  << inner_iter;
	break;
      }
      pdens_aux = ndens_aux;
    } while(dt_inner<dt_growth);
    dt_sgn /= 2.0;
    //std::cout << "\n[ee] Halving " << dt_sgn << " inner: " << dt_inner << " iter: " << inner_iter;
    if (inner_iter > max_iter) {
      std::cout << "\n[ee] Max iterations reached. Abort.\n" << inner_iter;
      std::terminate();
    }
    //std::cout << "\n[ee] Check while " << dt_inner << " < " << dt_growth << " iter: " << inner_iter;
  } while(/*(dt_inner<dt_growth)&&*/(negative_density));//while((dt_inner!=dt_growth) && (!negative_density));
  // update density
  ndens = ndens_aux;
  pdens = ndens;
  for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
    for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {	  
      srate2d[l][q] += srate2d_aux[l][q];
    }
  }
  //nrate2d += nrate2d_aux;
  return 1.0;
}

// ========================= SiH4 =========================

int NEvo::compute_split_sgnucleationsih4() {
  // Time splitting
  // + nucleation and growth at step dt/2
    
  double wsg = static_cast<double>(nm.rs.wsg);
  double wnu = static_cast<double>(nm.rs.wnu);

  unsigned int max_iter = 11;
  unsigned int inner_iter = 0;
  
  // delta time for particle growth (split)
  double dt_growth = 0.5 * nm.tm.ndeltat;
  
  // delta time for sg and nucleation
  double dt_sgn = dt_growth;
  double dt_inner = dt_growth;

  bool negative_density = false;

  int err;

  sih4rate = 0.0;
 
  do {
    inner_iter += 1;
    negative_density = false;
    pdens_aux = pdens;
    ndens_aux = zero2d;

    srate2d_aux = zero2d;
    nrate2d_aux = zero2d;
    dt_inner = 0.0;
    do {
      dt_inner += dt_sgn;
      //std::cout << "\n[ee] dt_sgn " << dt_sgn << " inner: " << dt_inner << " iter: " << inner_iter;

      sgrowth_total_rate = 0.0;
      if (wsg>0.0) {
	for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
	  for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
	    sgrowth_total_rate += (sgrowth_effective*vsih4
				   *pdens_aux[l][q]*surface_area[l]);
	  }
	}
      }
      nsih4 = nsih4_aux/(1.0+dt_sgn*(nm.rs.sih4nmol*nm.rs.nucleation_rate/nsih4_ini
				     + sgrowth_total_rate /*+ 10.981*/));//*
      sih4rate += (nsih4-nsih4_aux)/dt_sgn;

      err = compute_sgnuc_sih4();
      //std::cout << "\n[ee] Nucleation " <<  jnucleation[0][cr.gm.chrgs.maxnegative];
      for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
	for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
	  ndens_aux[l][q] = (pdens_aux[l][q] + dt_sgn
			  * (wsg*gsurfacegrowth[l][q]+wnu*jnucleation[l][q]));
	
	  // // store rates
	  // if (first_split) {
	  //   srate2d[l][q] = 0.5*wsg*gsurfacegrowth[l][q];
	  //   nrate2d[l][q] = 0.5*wnu*jnucleation[l][q];
	  // }
	  // else {
	  //   srate2d[l][q] += 0.5*wsg*gsurfacegrowth[l][q];
	  //   nrate2d[l][q] += 0.5*wnu*jnucleation[l][q];
	  // }
	  srate2d_aux[l][q] += 0.5*wsg*gsurfacegrowth[l][q];
	  //nrate2d_aux[l][q] += wnu*jnucleation[l][q];

	  if (ndens_aux[l][q]<0.0 && fabs(ndens_aux[l][q])<1.0) {
	    ndens_aux[l][q] = 0.0;
	  }
	  if (ndens_aux[l][q]<0.0){//EPSILON) {
	    std::cout << "\n[ee] Negative nanoparticle ndensity (growth + nuc) 1 ndens ("
		      << l << ", " << q << ") dt_sgn = " << dt_sgn
		      << " ndens[l][q] = " << ndens_aux[l][q];
	    //std::terminate();
	    negative_density = true;
	    dt_inner = 0.0;
	    break;
	  }
	  //gsurfacegrowth[l][q] = 0.0;
	}
	if (negative_density) {
	  std::cout << "\n[ee] Break 1 : " << inner_iter;
	  break;
	}
      }
      if (negative_density) {
	std::cout << "\n[ee] Break 2 : "  << inner_iter;
	break;
      }
      pdens_aux = ndens_aux;
      nsih4_aux = nsih4;
    } while(dt_inner<dt_growth);
    dt_sgn /= 2.0;
    //std::cout << "\n[ee] Halving " << dt_sgn << " inner: " << dt_inner << " iter: " << inner_iter;
    if (inner_iter > max_iter) {
      std::cout << "\n[ee] Max iterations reached. Abort.\n" << inner_iter;
      std::terminate();
    }
    //std::cout << "\n[ee] Check while " << dt_inner << " < " << dt_growth << " iter: " << inner_iter;
  } while(/*(dt_inner<dt_growth)&&*/(negative_density));//while((dt_inner!=dt_growth) && (!negative_density));
  // update density
  ndens = ndens_aux;
  pdens = ndens;

  // sgrowth_total_rate = 0.0;
  for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
    for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {	  
      srate2d[l][q] += srate2d_aux[l][q];
      // sgrowth_total_rate += srate2d[l][q];
    }
  }
  //nrate2d += nrate2d_aux;
  return 1.0;
}


// ========================= SiH4 =========================

int NEvo::compute_nosplit_sgnucleationsih4() {
    
  double wsg = static_cast<double>(nm.rs.wsg);
  double wnu = static_cast<double>(nm.rs.wnu);

  unsigned int max_iter = 11;
  unsigned int inner_iter = 0;
  
  int err;

  double dtg = nm.tm.ndeltat;
  
  sih4rate = 0.0;
 
  sgrowth_total_rate = 0.0;

  // WARNING no split
  pdens_aux = pdens;
  
  if (wsg>0.0) {
    for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
      for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
	sgrowth_total_rate += (sgrowth_effective*vsih4
			       *pdens[l][q]*surface_area[l]);
      }
    }
  }
  nsih4 = nsih4_aux/(1.0+dtg*(nm.rs.sih4nmol*nm.rs.nucleation_rate/nsih4_ini
				 + sgrowth_total_rate /*+ 10.981*/));//*
  sih4rate += (nsih4-nsih4_aux)/dtg;

  err = compute_sgnuc_sih4();

  nsih4_aux = nsih4;
  // sgrowth_total_rate = 0.0;
  for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
    for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {	  
      srate2d[l][q] = gsurfacegrowth[l][q]/dtg;
      // sgrowth_total_rate += srate2d[l][q];
    }
  }
  //nrate2d += nrate2d_aux;
  return 1.0;
}



// ------------------ openmp version ------------------

int NEvo::compute_sgrowth() {
//
// Explicit Surface Growth
//
// #pragma omp for collapse(2)
    for (unsigned int l = 1; l < cr.gm.vols.nsections-1; ++l) {
      for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
        gsurfacegrowth[l][q] = adim_srate[l-1]*pdens[l-1][q]
                             - adim_srate[l]*pdens[l][q];
        // Growth to last section only adds particles, removal requires to
        // add sections (note sign plus + )
        gsurfacegrowth[cr.gm.vols.nsections-1][q] = adim_srate[cr.gm.vols.nsections-2]
                                        * pdens[cr.gm.vols.nsections-2][q];
        // Growth of first section only removes particles
        // (note sign minus - )
        gsurfacegrowth[0][q] = -adim_srate[0]*pdens[0][q];
    }
  }
  return 0;
}

int NEvo::compute_sgrowth_adp() {
//
// Explicit Surface Growth
//
// #pragma omp for collapse(2)
    for (unsigned int l = 1; l < cr.gm.vols.nsections-1; ++l) {
      for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
        gsurfacegrowth[l][q] = adim_srate[l-1]*pdens_aux[l-1][q]
                             - adim_srate[l]*pdens_aux[l][q];
        // Growth to last section only adds particles, removal requires to
        // add sections (note sign plus + )
        gsurfacegrowth[cr.gm.vols.nsections-1][q] = adim_srate[cr.gm.vols.nsections-2]
                                        * pdens_aux[cr.gm.vols.nsections-2][q];
        // Growth of first section only removes particles
        // (note sign minus - )
        gsurfacegrowth[0][q] = -adim_srate[0]*pdens_aux[0][q];
    }
  }
  return 0;
}

// ========================= SiH4 =========================
int NEvo::compute_sih4() {
  
  return 0.0;
}

int NEvo::compute_sgnuc_sih4() {
//
// Explicit Surface Growth
//
  // for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
  //   for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
  //     jnucleation[l][q] = 0.0;
  //   }
  // }
  jnucleation[0][cr.gm.chrgs.maxnegative] = nm.rs.nucleation_rate*nsih4/nsih4_ini;


  // surface growth
  sgrowth_rate_sih4 = sgrowth_effective * nsih4 * vsih4 *sih4_vol;
  surface_rate = sgrowth_rate_sih4 * surface_area;
  
  adim_srate[0] = surface_rate[0] / (cr.gm.vols.volumes[1]-cr.gm.vols.volumes[0]);
  for (unsigned int l = 1; l < cr.gm.vols.nsections; ++l) {
    adim_srate[l] = surface_rate[l] / (cr.gm.vols.volumes[l+1]-cr.gm.vols.volumes[l]);
  }
  // Growth to last section only adds particles, removal requires to
  // add sections (note sign plus + )
  adim_srate[cr.gm.vols.nsections-1] = surface_rate[cr.gm.vols.nsections-2]
                           / (cr.gm.vols.volumes[cr.gm.vols.nsections-1]
                           -cr.gm.vols.volumes[cr.gm.vols.nsections-2]);
  
// #pragma omp for collapse(2)
  for (unsigned int l = 1; l < cr.gm.vols.nsections-1; ++l) {
    for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
      gsurfacegrowth[l][q] = adim_srate[l-1]*pdens_aux[l-1][q]
	- adim_srate[l]*pdens_aux[l][q];
      // Growth to last section only adds particles, removal requires to
      // add sections (note sign plus + )
      gsurfacegrowth[cr.gm.vols.nsections-1][q] = adim_srate[cr.gm.vols.nsections-2]
	* pdens_aux[cr.gm.vols.nsections-2][q];
      // Growth of first section only removes particles
      // (note sign minus - )
      gsurfacegrowth[0][q] = -adim_srate[0]*pdens_aux[0][q];
    }
  }
  return 0;
}

int NEvo::compute_coagulation() {

  //#pragma omp parallel num_threads(12)
  //{
  //#pragma omp parallel for schedule(auto) firstprivate(pdens)// Good
  for (unsigned int i=0; i<cr.eta_factor_vector.size(); ++i) {
    //
    EtaCreationFactor *ieta;
    ieta = &cr.eta_factor_vector[i];
//     #pragma omp atomic
    birth_vector[ieta->l_][ieta->q_] += ieta->eta_
                                      * pdens[ieta->m_][ieta->p_]
                                      * pdens[ieta->n_][ieta->r_];
  }

  // iterate in death_factor list
//   #pragma omp parallel for// schedule(auto)
//#pragma omp parallel for schedule(auto) firstprivate(pdens)// Good
  for (unsigned int j=0; j<cr.death_factor_vector.size(); ++j) {
    //
    DeathFactor *idth;
    idth = &cr.death_factor_vector[j];
//     #pragma omp atomic
    death_vector[idth->l_][idth->q_] += idth->death_
                                      * pdens[idth->m_][idth->p_]
                                      * pdens[idth->l_][idth->q_];
  }

  //}//parallel region
//   #pragma omp barrier
//   #pragma omp for collapse(2)
  for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
    for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
      kcoagulation[l][q] = birth_vector[l][q]
                         - death_vector[l][q];

      // WARNING DANGER
//       if(abs(kcoagulation[l][q])<EPSILONETA) {
//         kcoagulation[l][q] = 0.0;
// //         std::cout << "\n[ee] Negative kcoagulation";
// //         std::terminate();
//       }
      // set to zero
      birth_vector[l][q]=0.0; death_vector[l][q]=0.0;
    }
  }
  
  return 0;
}

// ------------------ openmp version ------------------

int NEvo::compute_coagulation_omp() {
  // BOOST_LOG_SEV(lg, info) << "COMP";
  double bsum = 0.0;
  double dsum = 0.0;
  for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
    for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
      birth_vector[l][q]=0.0; death_vector[l][q]=0.0;
    }
  }
#pragma omp parallel
  {
#pragma omp for nowait// schedule(dynamic, 100)// schedule(static, 1)// schedule(guided)//, 32768)// firstprivate(pdens)// Good
    for (unsigned int i=0; i<cr.eta_factor_vector.size(); ++i) {
      unsigned int id;
      short l;
      short q;
      short m;
      short p;
      short n;
      short r;
      double eta;
      //get_eta(cr.eta_factor_vector[i], id, l, q, m, p, n, r, eta);
      cr.eta_factor_vector[i].get_eta2(id, l, q, m, p, n, r, eta);
      birth_vector[l][q] += eta * (pdens[m][p] * pdens[n][r]);
      //beans[i] += eta * (pdens[m][p] * pdens[n][r]);
      //birth_vector[l][q] = 0.0;
    }

    // iterate in death_factor list
#pragma omp for nowait// schedule(dynamic, 100)// schedule(static, 1)// schedule(guided)//, 32768)// firstprivate(pdens)// Good
    for (unsigned int j=0; j<cr.death_factor_vector.size(); ++j) {     
      unsigned int id;
      short l;
      short q;
      short m;
      short p;
      double death;
      //get_death(cr.death_factor_vector[i], id, l, q, m, p, death);
      cr.death_factor_vector[j].get_death2(id, l, q, m, p, death);
      death_vector[l][q] += death * pdens[m][p];// (pdens[m][p] * pdens[l][q]);
      //death_vector[l][q] = 0.0;
      //beans += death * (pdens[m][p] * pdens[l][q]);
    }
//     // iterate in death_factor list
// #pragma omp for nowait// schedule(static, 1)// schedule(guided)//, 32768)// firstprivate(pdens)// Good
//     for (unsigned int j=0; j<cr.death_factor_vector.size(); ++j) {     
//       unsigned int id;
//       short l;
//       short q;
//       short m;
//       short p;
//       double death;
//       //get_death(cr.death_factor_vector[i], id, l, q, m, p, death);
//       cr.death_factor_vector[j].get_death2(id, l, q, m, p, death);
//       death_vector[l][q] += death * (pdens[m][p] * pdens[l][q]);
//       unsigned int i = j;
//       if(i < cr.eta_factor_vector.size()) {
// 	short n;
// 	short r;
// 	double eta;
// 	cr.eta_factor_vector[i].get_eta2(id, l, q, m, p, n, r, eta);
// 	birth_vector[l][q] += eta * (pdens[m][p] * pdens[n][r]);
//       }
//       //death_vector[l][q] = 0.0;
//       //beans += death * (pdens[m][p] * pdens[l][q]);
//     }


#pragma omp barrier
  
#pragma omp for collapse(2)
    for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
      for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
	kcoagulation[l][q] = birth_vector[l][q]
	  - death_vector[l][q]*pdens[l][q];

	// WARNING DANGER
	//       if(abs(kcoagulation[l][q])<EPSILONETA) {
	//         kcoagulation[l][q] = 0.0;
	// //         std::cout << "\n[ee] Negative kcoagulation";
	// //         std::terminate();
	//       }
	// set to zero
	//birth_vector[l][q]=0.0; death_vector[l][q]=0.0;
      }
    }
  
    
// #pragma omp barrier

// #pragma omp for collapse(2) reduction(+:bsum,dsum)
//     for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
//       for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
// 	bsum += birth_vector[l][q];
// 	dsum += death_vector[l][q];
// 	birth_vector[l][q]=0.0; death_vector[l][q]=0.0;
//       }
//     }

// #pragma omp single
//     {
//       if(is_close(bsum-dsum, 0.0)==false) {
// 	BOOST_LOG_SEV(lg, info) << "Volume not conserved";
// 	BOOST_LOG_SEV(lg, info) << bsum;
// 	BOOST_LOG_SEV(lg, info) << dsum;
// 	BOOST_LOG_SEV(lg, info) << dsum/bsum;
// 	std::terminate();
//       }
//       else {
// 	BOOST_LOG_SEV(lg, info) << "birth and death";
// 	BOOST_LOG_SEV(lg, info) << bsum;
// 	BOOST_LOG_SEV(lg, info) << dsum;
//       }
//     }
  }
  return 0;
}
