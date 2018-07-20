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

#include "CRate.h"

int CRate::open() {
  std::string fullname = dirname + prefix_filename + ".h5";

  int err = open_hdf5(fullname, h5obj, "write");
  if(err!=0) {
    BOOST_LOG_SEV(lg, error) << "h5 I/O error opening file. Terminate.";
    std::terminate();
  }
  BOOST_LOG_SEV(lg, info) << "Open file for read/write: "
                          << fullname << " -> Success. Id: "
                          << h5obj.getId();
  return 0;
}

int CRate::close() {
  hid_t id = h5obj.getId();

  int err = close_hdf5(h5obj);
  if(err!=0) {
    BOOST_LOG_SEV(lg, error) << "h5 I/O error closing file. Terminate.";
    std::terminate();
  }
  BOOST_LOG_SEV(lg, info) << "Closed file Id: " << id;
  return 0;
}

int CRate::read() {
  hid_t id = h5obj.getId();

  gm = GridModel(h5obj, lg);

  gm.read();

  // grid definition
  grid = {{gm.vols.nsections, gm.chrgs.nsections}};
  grid4 = {{gm.vols.nsections, gm.chrgs.nsections,
            gm.vols.nsections, gm.chrgs.nsections}};

  return 0;
}

int CRate::read_results() {

  int err;

  BOOST_LOG_SEV(lg, info) << "Reading grid";
  err = read();

  BOOST_LOG_SEV(lg, info) << "Reading eta creation factor";
  read_etafactor();

  BOOST_LOG_SEV(lg, info) << "Reading death factor";
  read_deathfactor();

  return 0;
}

int CRate::write() {
  hid_t id = h5obj.getId();

  BOOST_LOG_SEV(lg, info) << "Writing eta creation factor";
  write_etafactor();

  BOOST_LOG_SEV(lg, info) << "Writing death factor";
  write_deathfactor();

  BOOST_LOG_SEV(lg, info) << "Writing potentials";
  write_potentials();
  
  BOOST_LOG_SEV(lg, info) << "Writing enhancement factor";
  write_efactor();

  BOOST_LOG_SEV(lg, info) << "Writing coagulation rate";
  write_rcoagulation();

  return 0;
}

int CRate::write_frompairs() {
  hid_t id = h5obj.getId();

  BOOST_LOG_SEV(lg, info) << "Writing eta creation factor";
  write_etafactor();

  BOOST_LOG_SEV(lg, info) << "Writing death factor";
  write_deathfactor();

  BOOST_LOG_SEV(lg, info) << "Writing potentials";
  write_potentials();
  
  BOOST_LOG_SEV(lg, info) << "Writing enhancement factor!!";
  write_efactor_serial();

  BOOST_LOG_SEV(lg, info) << "Writing potentials";
  write_potentials_serial();
  
  BOOST_LOG_SEV(lg, info) << "Writing coagulation rate!!";
  write_rcoagulation();

  return 0;
}



int CRate::compute_list(std::vector<std::string> sfilelist_) {
  //
  BOOST_LOG_SEV(lg, info) << "Processing pairs files...";
  auto start = std::chrono::system_clock::now();

  sfilelist = sfilelist_;

  BOOST_LOG_SEV(lg, info) << "Reading h5 list";
  
  for(unsigned int isf=0; isf<sfilelist.size(); ++isf) {

    H5::H5File h5fl;
     
    std::string fullname = dirname + sfilelist[isf];
    
    int err = open_hdf5(fullname, h5fl, "write");
    if(err!=0) {
      BOOST_LOG_SEV(lg, error) << "h5 I/O error opening file. Terminate.";
      std::terminate();
    }
    hid_t id = h5fl.getId();
    BOOST_LOG_SEV(lg, info) << "Open file for read/write: "
			    << fullname << " -> Success. Id: "
			    << id;

    std::string gname = "Enhancement_factor_serial";
    std::string dsname = "Indices";
    boost_short_array2d efindices;
    read_dset2d_hdf5<boost_short_array2d, short>(h5fl, gname, dsname, efindices);

    std::string edsname = "efactor";
    darray daefactor;
    read_dset_hdf5<darray, double>(h5fl, gname, edsname, daefactor);

    // grid definition
    grid = {{gm.vols.nsections, gm.chrgs.nsections}};
    grid4 = {{gm.vols.nsections, gm.chrgs.nsections,
	      gm.vols.nsections, gm.chrgs.nsections}};
  
    // resizing
    efactor.resize(grid4);
    cpotentials.resize(grid4);
    bpotentials.resize(grid4);
    rbarriers.resize(grid4);
    rcoag.resize(grid4);
   

    for (unsigned int i = 0; i<efindices.shape()[0]; ++i) {
      short l = efindices[i][0];
      short q = efindices[i][1];
      short m = efindices[i][2];
      short p = efindices[i][3];
      double etaf = daefactor[i];
      efactor[l][q][m][p] = etaf;
    }
    /////////

    std::string cpgname = "Contact_potential_serial";
    std::string cpidsname = "Indices";
    boost_short_array2d cpindices;
    read_dset2d_hdf5<boost_short_array2d, short>(h5fl, cpgname, cpidsname, cpindices);

    std::string cpdsname = "contact_potential";
    darray dacpotentials;
    read_dset_hdf5<darray, double>(h5fl, cpgname, cpdsname, dacpotentials);
   

    for (unsigned int i = 0; i<cpindices.shape()[0]; ++i) {
      short l = cpindices[i][0];
      short q = cpindices[i][1];
      short m = cpindices[i][2];
      short p = cpindices[i][3];
      double cpot = dacpotentials[i];
      cpotentials[l][q][m][p] = cpot;
    }
    ////////

    std::string bpgname = "Barrier_potential_serial";
    std::string bpidsname = "Indices";
    boost_short_array2d bpindices;
    read_dset2d_hdf5<boost_short_array2d, short>(h5fl, bpgname, bpidsname, bpindices);

    std::string bpdsname = "barrier_potential";
    darray dabpotentials;
    read_dset_hdf5<darray, double>(h5fl, bpgname, bpdsname, dabpotentials);

    std::string bcpdsname = "contact_potential";
    darray dabcpotentials;
    read_dset_hdf5<darray, double>(h5fl, bpgname, bcpdsname, dabcpotentials);

    std::string rbdsname = "rbarrier";
    darray darbarriers;
    read_dset_hdf5<darray, double>(h5fl, bpgname, rbdsname, darbarriers);

    for (unsigned int i = 0; i<bpindices.shape()[0]; ++i) {
      short l = bpindices[i][0];
      short q = bpindices[i][1];
      short m = bpindices[i][2];
      short p = bpindices[i][3];
      double bpot = dabpotentials[i];
      bpotentials[l][q][m][p] = bpot;
      double cpot = dabcpotentials[i];
      cpotentials[l][q][m][p] = cpot;
      double rb = darbarriers[i];
      rbarriers[l][q][m][p] = rb;
    }        

    ///////
    err = close_hdf5(h5fl);
    BOOST_LOG_SEV(lg, info) << "Closed file Id: " << id;
    
  }

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double>  elapsed_seconds = end-start;
  BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();

  // Compute beta0
  beta0 = pow(3.0/(4.0*M_PI), 1.0/6.0)
        * sqrt(6.0*Kboltz*gm.gsys.temperature/gm.gsys.nmdensity);

  start = std::chrono::system_clock::now();
  BOOST_LOG_SEV(lg, info) << "Computing coagulation rate";
  compute_rcoagulation();
  BOOST_LOG_SEV(lg, info) << "Done... coagulation rate";
  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();
  return 0;
}

void CRate::bind_efactorfunc(){
  // bind efactorfunc
  switch(gm.einter.method) {
    case 0:
      efactorfunc = &CRate::compute_efactor_mpc;
      BOOST_LOG_SEV(lg, info) << "Bounded method MPC : " << gm.einter.method;
      break;
    case 1:
      efactorfunc = &CRate::compute_efactor_ipa;
      BOOST_LOG_SEV(lg, info) << "Bounded method IPA : " << gm.einter.method;
      break; 
    case 2:
      efactorfunc = &CRate::compute_efactor_coul;
      BOOST_LOG_SEV(lg, info) << "Bounded method Coulomb : " << gm.einter.method;
      break;
    default:
      BOOST_LOG_SEV(lg, error) << "Method " << gm.einter.method
                               << " not known";
      std::terminate();
  }
}


int CRate::write_pairs() {

  BOOST_LOG_SEV(lg, info) << "Computing and writing particle pairs...";
  
  // grid definition
  grid = {{gm.vols.nsections, gm.chrgs.nsections}};
  grid4 = {{gm.vols.nsections, gm.chrgs.nsections,
            gm.vols.nsections, gm.chrgs.nsections}};
  
  // resizing
  efactor.resize(grid4);
  cpotentials.resize(grid4);
  bpotentials.resize(grid4);
  rbarriers.resize(grid4);
  rcoag.resize(grid4);

  BOOST_LOG_SEV(lg, info) << "Size for array efactor : " << efactor.size();
  // BOOST_LOG_SEV(lg, info) <<"\n r size " << gm.vols.radii.size();
  // std::cerr << "\n q size " << gm.chrgs.charges.size();
  BOOST_LOG_SEV(lg, info) << "Radii size : " << gm.vols.nsections;
  BOOST_LOG_SEV(lg, info) << "Charge size : " << gm.chrgs.nsections;
  
  enhancement::Enhancement enh(gm.vols.radii,
			       gm.chrgs.charges*eCharge,
			       efactor,
			       cpotentials,
			       bpotentials,
			       rbarriers,
			       gm.einter.dconstant,
			       lg);

  BOOST_LOG_SEV(lg, info) << "Computing reduced particle pairs...";
  auto start = std::chrono::system_clock::now();
  //
  enh.compute_reducedpairs();
  //
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();


  std::string fullname = dirname + prefix_filename;
  
  BOOST_LOG_SEV(lg, info) << "Writing particle pairs...";
  start = std::chrono::system_clock::now();
  //
  enh.write_particlepairs(fullname);
  //
  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();

  return 0;
}


int CRate::read_pairs() {

  BOOST_LOG_SEV(lg, info) << "Reading particle pairs...";
  
  // grid definition
  grid = {{gm.vols.nsections, gm.chrgs.nsections}};
  grid4 = {{gm.vols.nsections, gm.chrgs.nsections,
            gm.vols.nsections, gm.chrgs.nsections}};
  
  // resizing
  efactor.resize(grid4);
  cpotentials.resize(grid4);
  bpotentials.resize(grid4);
  rbarriers.resize(grid4);
  rcoag.resize(grid4);

  BOOST_LOG_SEV(lg, info) << "Size for array efactor : " << efactor.size();
  // BOOST_LOG_SEV(lg, info) <<"\n r size " << gm.vols.radii.size();
  // std::cerr << "\n q size " << gm.chrgs.charges.size();
  BOOST_LOG_SEV(lg, info) << "Radii size : " << gm.vols.nsections;
  BOOST_LOG_SEV(lg, info) << "Charge size : " << gm.chrgs.nsections;
  
  enhancement::Enhancement enh(gm.vols.radii,
			       gm.chrgs.charges*eCharge,
			       efactor,
			       cpotentials,
			       bpotentials,
			       rbarriers,
			       gm.einter.dconstant,
			       lg);

  std::string fullname = dirname + prefix_filename;
  
  BOOST_LOG_SEV(lg, info) << "Reading particle pairs...";
  auto start = std::chrono::system_clock::now();
  //
  enh.read_particlepairs(fullname);
  //
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double>  elapsed_seconds = end-start;
  BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();

  return 0;
}

int CRate::compute_frompairs() {

  BOOST_LOG_SEV(lg, info) << "Reading particle pairs...";
  
  BOOST_LOG_SEV(lg, info) << "Radii size : " << gm.vols.nsections;
  BOOST_LOG_SEV(lg, info) << "Charge size : " << gm.chrgs.nsections;
  
    
  enhancement::Enhancement enh(gm.vols.radii,
			       gm.chrgs.charges*eCharge,
  			       gm.einter.dconstant,
			       lg);

  std::string fullname = dirname + prefix_filename;
  
  BOOST_LOG_SEV(lg, info) << "Reading particle pairs...";
  auto start = std::chrono::system_clock::now();
  //
  enh.read_particlepairs(fullname);
  //
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double>  elapsed_seconds = end-start;
  BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();


  switch(gm.einter.method) {
    case 0:
      BOOST_LOG_SEV(lg, info) << "Bounded method MPC : " << gm.einter.method;
      BOOST_LOG_SEV(lg, info) << "Computing MPC potentials at contact...";
      start = std::chrono::system_clock::now();
      //
      enh.compute_mpcpotential_contact();
      //
      end = std::chrono::system_clock::now();
      elapsed_seconds = end-start;
      BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();
      BOOST_LOG_SEV(lg, info) << "Computing MPC potentials barriers...";
      start = std::chrono::system_clock::now();
      //
      enh.compute_mpcpotential_barrier();
      //
      end = std::chrono::system_clock::now();
      elapsed_seconds = end-start;
      BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();
      break;
    case 1:
      BOOST_LOG_SEV(lg, info) << "Bounded method IPA : " << gm.einter.method;
      BOOST_LOG_SEV(lg, info) << "Computing IPA potentials at contact...";
      start = std::chrono::system_clock::now();
      //
      enh.compute_ipapotential_contact();
      //
      end = std::chrono::system_clock::now();
      elapsed_seconds = end-start;
      BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();
      BOOST_LOG_SEV(lg, info) << "Computing IPA potentials barriers...";
      start = std::chrono::system_clock::now();
      //      
      enh.compute_ipapotential_barrier();
      //
      end = std::chrono::system_clock::now();
      elapsed_seconds = end-start;
      BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();
      break; 
    case 2:
      BOOST_LOG_SEV(lg, info) << "Bounded method Coulomb : " << gm.einter.method;
      BOOST_LOG_SEV(lg, info) << "Computing coulomb potentials at contact...";
      start = std::chrono::system_clock::now();
      //      
      enh.compute_coulombpotential_contact();
      //
      end = std::chrono::system_clock::now();
      elapsed_seconds = end-start;
      BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();
      break;
  default:
    BOOST_LOG_SEV(lg, error) << "Method " << gm.einter.method
  			     << " not known";
      std::terminate();
  }
 
  BOOST_LOG_SEV(lg, info) << "Computing enhancement factor grid";
  start = std::chrono::system_clock::now();
  //
  enh.compute_enhancementfactor_frompairs();
  //
  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();

  enh.get_efindices(daefactor, efindices);

  enh.get_cpindices(dacpotentials, cpindices);
    
  enh.get_bpindices(dabcpotentials, dabpotentials, darbarriers, bpindices);
  //for (unsigned int ii=0; ii<daefactor.size(); ++ii) {
  //  std::cout << '\n' << daefactor[ii];
  //}
  
  // Compute beta0
  beta0 = pow(3.0/(4.0*M_PI), 1.0/6.0)
        * sqrt(6.0*Kboltz*gm.gsys.temperature/gm.gsys.nmdensity);

  // // BOOST_LOG_SEV(lg, info) << "Computing enhancement factor grid";
  // // compute_efactor_grid();
  // // BOOST_LOG_SEV(lg, info) << "Done... ehancement factor grid";

  // for a large grid 100x300 comment lines below
  start = std::chrono::system_clock::now();
  BOOST_LOG_SEV(lg, info) << "Computing coagulation rate";
  compute_rcoagulation();
  BOOST_LOG_SEV(lg, info) << "Done... coagulation rate";
  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();
  

  BOOST_LOG_SEV(lg, info) << "Computing eta creation rate";
  compute_etafactor();
  BOOST_LOG_SEV(lg, info) << "Done... eta";

  BOOST_LOG_SEV(lg, info) << "Computing death rate";
  compute_deathfactor();
  BOOST_LOG_SEV(lg, info) << "Done... death";

  return 0;
}


int CRate::compute() {

  // grid definition
  grid = {{gm.vols.nsections, gm.chrgs.nsections}};
  grid4 = {{gm.vols.nsections, gm.chrgs.nsections,
            gm.vols.nsections, gm.chrgs.nsections}};
  
  // resizing
  efactor.resize(grid4);
  cpotentials.resize(grid4);
  bpotentials.resize(grid4);
  rbarriers.resize(grid4);
  rcoag.resize(grid4);

  BOOST_LOG_SEV(lg, info) << "Size for array efactor : " << efactor.size();
  BOOST_LOG_SEV(lg, info) << "Radii size : " << gm.vols.nsections;
  BOOST_LOG_SEV(lg, info) << "Charge size : " << gm.chrgs.nsections;     

  enhancement::Enhancement enh(gm.vols.radii,
			       gm.chrgs.charges*eCharge,
			       efactor,
			       cpotentials,
			       bpotentials,
			       rbarriers,
			       gm.einter.dconstant,
			       lg);

  BOOST_LOG_SEV(lg, info) << "Computing reduced particle pairs...";
  auto start = std::chrono::system_clock::now();
  //
  enh.compute_reducedpairs();
  //
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();  
  
  switch(gm.einter.method) {
    case 0:
      BOOST_LOG_SEV(lg, info) << "Bounded method MPC : " << gm.einter.method;
      BOOST_LOG_SEV(lg, info) << "Computing MPC potentials at contact...";
      start = std::chrono::system_clock::now();
      //
      enh.compute_mpcpotential_contact();
      //
      end = std::chrono::system_clock::now();
      elapsed_seconds = end-start;
      BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();
      BOOST_LOG_SEV(lg, info) << "Computing MPC potentials barriers...";
      start = std::chrono::system_clock::now();
      //
      enh.compute_mpcpotential_barrier();
      //
      end = std::chrono::system_clock::now();
      elapsed_seconds = end-start;
      BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();
      break;
    case 1:
      BOOST_LOG_SEV(lg, info) << "Bounded method IPA : " << gm.einter.method;
      BOOST_LOG_SEV(lg, info) << "Computing IPA potentials at contact...";
      start = std::chrono::system_clock::now();
      //
      enh.compute_ipapotential_contact();
      //
      end = std::chrono::system_clock::now();
      elapsed_seconds = end-start;
      BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();
      BOOST_LOG_SEV(lg, info) << "Computing IPA potentials barriers...";
      start = std::chrono::system_clock::now();
      //      
      enh.compute_ipapotential_barrier();
      //
      end = std::chrono::system_clock::now();
      elapsed_seconds = end-start;
      BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();
      break; 
    case 2:
      BOOST_LOG_SEV(lg, info) << "Bounded method Coulomb : " << gm.einter.method;
      BOOST_LOG_SEV(lg, info) << "Computing coulomb potentials at contact...";
      start = std::chrono::system_clock::now();
      //      
      enh.compute_coulombpotential_contact();
      //
      end = std::chrono::system_clock::now();
      elapsed_seconds = end-start;
      BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();
      break;
    case 3:
      BOOST_LOG_SEV(lg, info) << "Bounded Hybrid method : " << gm.einter.method;
      BOOST_LOG_SEV(lg, info) << "Computing MPC potentials at contact...";
      start = std::chrono::system_clock::now();
      //
      enh.compute_mpcpotential_contact();
      //
      end = std::chrono::system_clock::now();
      elapsed_seconds = end-start;
      BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();
      BOOST_LOG_SEV(lg, info) << "Computing IPA potentials barriers...";
      start = std::chrono::system_clock::now();
      //
      enh.compute_ipapotential_barrier();
      //
      end = std::chrono::system_clock::now();
      elapsed_seconds = end-start;
      BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();
      break;      
  default:
    BOOST_LOG_SEV(lg, error) << "Method " << gm.einter.method
  			     << " not known";
      std::terminate();
  }
 
  BOOST_LOG_SEV(lg, info) << "Computing enhancement factor grid";
  start = std::chrono::system_clock::now();
  //
  enh.compute_enhancement_factor();
  //
  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();
  
  // // Compute electrostatic to thermal ratio
  // electhermratio = 1.0*eCharge*eCharge
  //                 / (4.0*M_PI*EpsilonZero*Kboltz*gm.gsys.temperature);

  // BOOST_LOG_SEV(lg, info) << "Electrostatic to thermal ratio: "
  //                         << electhermratio;

  // Compute beta0
  beta0 = pow(3.0/(4.0*M_PI), 1.0/6.0)
        * sqrt(6.0*Kboltz*gm.gsys.temperature/gm.gsys.nmdensity);

  //BOOST_LOG_SEV(lg, info) << "Computing enhancement factor grid";
  //compute_efactor_grid();
  //BOOST_LOG_SEV(lg, info) << "Done... ehancement factor grid";

  start = std::chrono::system_clock::now();
  BOOST_LOG_SEV(lg, info) << "Computing coagulation rate";
  compute_rcoagulation();
  BOOST_LOG_SEV(lg, info) << "Done... coagulation rate";
  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();
  
  BOOST_LOG_SEV(lg, info) << "Computing eta creation rate";
  compute_etafactor();
  BOOST_LOG_SEV(lg, info) << "Done... eta";

  BOOST_LOG_SEV(lg, info) << "Computing death rate";
  compute_deathfactor();
  BOOST_LOG_SEV(lg, info) << "Done... death";

  return 0;
}

int CRate::compute_sym() {

  // bind efactorfunc
  bind_efactorfunc();

  // grid definition
  grid = {{gm.vols.nsections, gm.chrgs.nsections}};
  grid4 = {{gm.vols.nsections, gm.chrgs.nsections,
            gm.vols.nsections, gm.chrgs.nsections}};

  // resizing
  efactor.resize(grid4); 
  rcoag.resize(grid4);

  // Compute electrostatic to thermal ratio
  electhermratio = 1.0*eCharge*eCharge
                  / (4.0*M_PI*EpsilonZero*Kboltz*gm.gsys.temperature);

  BOOST_LOG_SEV(lg, info) << "Electrostatic to thermal ratio: "
                          << electhermratio;

  // Compute beta0
  beta0 = pow(3.0/(4.0*M_PI), 1.0/6.0)
        * sqrt(6.0*Kboltz*gm.gsys.temperature/gm.gsys.nmdensity);

  BOOST_LOG_SEV(lg, info) << "Computing enhancement factor grid";
  compute_efactor_grid_sym();
  BOOST_LOG_SEV(lg, info) << "Done... ehancement factor grid";

  BOOST_LOG_SEV(lg, info) << "Computing coagulation rate";
  compute_rcoagulation();
  BOOST_LOG_SEV(lg, info) << "Done... coagulation rate";

  BOOST_LOG_SEV(lg, info) << "Computing eta creation rate";
  compute_etafactor();
  BOOST_LOG_SEV(lg, info) << "Done... eta";

  BOOST_LOG_SEV(lg, info) << "Computing death rate";
  compute_deathfactor();
  BOOST_LOG_SEV(lg, info) << "Done... death";
  
  return 0;
}


void CRate::compute_efactor_grid_sym() {
  BOOST_LOG_SEV(lg, info) << "Computing symmetric version";
  // iterate in non repeated combinations of particle pairs (r,q)
#pragma omp parallel for collapse(4) schedule(auto)
  for (unsigned int l=0; l<gm.vols.nsections; ++l) {
    // iterate in charges particle 1
    for (unsigned int q=0; q<gm.chrgs.nsections; ++q) {
      // iterate in radii particle 2
      for (unsigned int m=0; m<gm.vols.nsections; ++m) {
	// iterate in charges particle 2	
	for (unsigned int p=0; p<gm.chrgs.nsections; ++p) {
	  unsigned int p1 = q * gm.vols.nsections + l;
	  unsigned int p2 = p * gm.vols.nsections + m;
	  // avoid repetitions
	  if((p>=q) && (p2 >= p1)){
	    efactor[l][q][m][p] = compute_efactor(gm.vols.radii[l],
						  gm.vols.radii[m],
						  gm.chrgs.charges[q],
						  gm.chrgs.charges[p]);
	    // write symmetric enhancement factor
	    efactor[m][p][l][q] =  efactor[l][q][m][p]; 
	  }
	}
      }
    }
  }
}

void CRate::compute_efactor_grid() {
  // iterate in radii particle 1
  //#pragma omp parallel for default(none) collapse(4) schedule(auto)
#pragma omp parallel for collapse(4) schedule(auto)
  for (unsigned int l=0; l<gm.vols.nsections; ++l) {
    // iterate in charges particle 1
    for (unsigned int q=0; q<gm.chrgs.nsections; ++q) {
      // iterate in radii particle 2
      for (unsigned int m=0; m<gm.vols.nsections; ++m) {
        // iterate in charges particle 2	
        for (unsigned int p=0; p<gm.chrgs.nsections; ++p) {
          efactor[l][q][m][p] = compute_efactor(gm.vols.radii[l],
                                                gm.vols.radii[m],
                                                gm.chrgs.charges[q],
                                                gm.chrgs.charges[p]);
        }
      }
    }
  }
}

inline 
double CRate::beta_free(const double vol1,
                        const double vol2,
                        const double beta0) {
  return beta0 * std::sqrt(1.0/vol1 +1.0/vol2)
         * std::pow(std::pow(vol1, 1.0/3.0)+std::pow(vol2, 1.0/3.0), 2);
}

void CRate::compute_rcoagulation() {
// WARNING beta symmetric optimization FIXME
// #pragma omp parallel for schedule(auto)
  for (unsigned int l = 0; l < gm.vols.nsections; ++l) {
    for (unsigned int q = 0; q < gm.chrgs.nsections; ++q) {
      for (unsigned int m = 0; m < gm.vols.nsections; ++m) {
        for (unsigned int p = 0; p < gm.chrgs.nsections; ++p) {
          rcoag[l][q][m][p] = gm.einter.multiplier * efactor[l][q][m][p]
                            * beta_free(gm.vols.volumes[l],
                                        gm.vols.volumes[m],
                                        beta0);

          if(rcoag[l][q][m][p]<0.0) {
            BOOST_LOG_SEV(lg, error) << "Negative coagulation rate. Terminate";
            std::terminate();
          }
        }
      }
    }
  }

}

// only volume pivoting
void CRate::compute_etafactor() {
  // eta function
  double eta;
  double comp1_sum;
  double charge_sum;
//   double max_charge1; // max charge for sum volume vm + vn
//   double max_charge2; // max charge for sum volume vm + vn
  unsigned int isum = 0;

  double eta_fraction;

  EtaCreationFactor eta_factor;

  // section boundaries for volumes
  darray vifaces = gm.vols.interfaces;

  // volume pivots
  darray vpivots = gm.vols.volumes;

  darray qpivots = gm.chrgs.charges;


  double lower_comp1, upper_comp1, curr_comp1;
  double curr_charge;

  isum = 0;
  // loops in l from 1st section
// #pragma omp parallel for schedule(auto) firstprivate(vifaces, vpivots, qpivots, rcoag) doesnt work
  for (unsigned int l = 0; l < gm.vols.nsections; ++l) {
    // HACK
    // special case for first section l=0
    //     p0    p   p1
    // |....*....|....*....|
    // i0        i1        i2
    // v in Zone  II  (v-v0)/(v1-v0)
    // in this case there is not lower volume
    if (l == 0) {
      // artificially rules out the comparison in the nonexistent zone I
      // v2 <= v <= v1
      lower_comp1 = vifaces[0];

    }
    else {
      lower_comp1 = vpivots[l-1];
    }

    // pivot at center
    curr_comp1 = vpivots[l];

    // special case for last section i=vols.nsecs1-1 sN-1
    //
    //           pN-2      p      pN-1                Pivots
    // ....|.......*.......|........*........|
    //    iN-2            iN-1              iN      Interfaces
    //
    // p must be in Zone I (pN-1 - p)/(pN-1 - pN-2)
    // In this case there is not lower bound
    if (l == gm.vols.nsections-1) {
      // artificially rules out the comparison in the nonexistent zone II
      // vN-1 <= v <= vN -> vN-1 <= v <= vN-2
      upper_comp1 = vifaces[l+1];
    }
    else {
      upper_comp1 = vpivots[l+1];
    }
    for (unsigned int q = 0; q < gm.chrgs.nsections; ++q) {
      curr_charge = qpivots[q];
      // particle with comp1 m and curr_charge p coalesce with comp1 n and curr_charge r
      // (m, p) + (n, r) -> (l, q) then comp1=m+n and charge=r+p
      // loops in volume sections m[0, l]
      for (unsigned int m = 0; m < l+1; ++m) {
//       for (unsigned int m = 0; m < vols.nsecs; ++m) {
        // loops in charges p[0, q]
        for (unsigned int p = 0; p < gm.chrgs.nsections; ++p) {
// DANGER
          for (unsigned int n = 0; n < l+1; ++n) {
//           for (unsigned int n = 0; n < vols.nsecs; ++n) {
//           for (unsigned int n = m; n < vols.nsecs; ++n) {
            // loops in charges
            for (unsigned int r = 0; r < gm.chrgs.nsections; ++r) {
//             for (unsigned int r = qind[n][1]; r > qind[n][0]; --r) {
            comp1_sum = vpivots[m] + vpivots[n];
            charge_sum = qpivots[r] + qpivots[p];
//             std::cerr << "\n[WW] sum , piv : " << charge_sum << '\t' << qifaces[q];
            if(charge_sum == curr_charge) {// WARNING float comparison FIXME below
//             if(logically_equal(charge_sum, curr_charge, 0.1)) {
//             if((abs(comp2_sum - curr_comp2) <= 0.01)) {
//               std::cerr << "\n[EE] equals p,r : " << p << '\t' << r;
              // vl-1 <= v < vl and Qq-1 <= Q < Qq

              if ((comp1_sum >= lower_comp1) && (comp1_sum < curr_comp1)) {
                //
                eta_fraction = (comp1_sum-lower_comp1)
                             / (curr_comp1-lower_comp1);
                // TODO abs not needed nor dirac p r
//                 eta = abs(eta_fraction) * (1.0 - 0.5 * kron_delta(m, n)
//                    /* * kron_delta(p,r)*/)*rcoag[l][q][m][p];
                // WARNING abs unnecessary eta_fraction > 0
//                 eta = (1.0-0.5*kron_delta(m, n))*/*0.5**/abs(eta_fraction) * rcoag[l][q][m][p];
                eta = 0.5 * eta_fraction * rcoag[l][q][m][p];
                if(eta > EPSILONETA) {
                  fill_eta(eta_factor, isum, l, q, m, p, n, r, eta);
                  eta_factor_vector.push_back(eta_factor);
                  isum++;
                } else {
                  if(eta<0.0) std::terminate;
                }
              }
              // vl <= v < vl+1 and Qq <= Q < Qq+1
              if ((comp1_sum >= curr_comp1) && (comp1_sum < upper_comp1)) {
                eta_fraction = (upper_comp1-comp1_sum)
                             / (upper_comp1-curr_comp1);
                // TODO abs not needed nor dirac p r
//                 eta = abs(eta_fraction) * (1.0 - 0.5 * kron_delta(m, n)
//                     /** kron_delta(p,r)*/)*rcoag[l][q][m][p];
//                 eta = (1.0-0.5*kron_delta(m, n))*/* 0.5**/abs(eta_fraction) * rcoag[l][q][m][p];
                eta = 0.5 * eta_fraction * rcoag[l][q][m][p];//Kumar pag. 221
                // WARNING EPSILONETA
                if(eta > EPSILONETA) {
                  fill_eta(eta_factor, isum, l, q, m, p, n, r, eta);
                  eta_factor_vector.push_back(eta_factor);
                  isum++;
                } else {
                  if(eta<0.0) std::terminate;
                }
              }
            } // end charge comp
           } // end for n
          } // end for r
        } // end for p
      } // end for m
    } // end for q charges
  } //end for l volumes

  BOOST_LOG_SEV(lg, info) << "eta nonzeros = " << isum;

  double max_eta = eta_factor_vector[0].eta_;
  double min_eta = eta_factor_vector[0].eta_;
  double avg_eta = 0.0;
  // iterate in eta_factor list
//   #pragma omp parallel for shared(max_eta, min_eta)// schedule(auto)
//   std::cout << std::endl << "Eta " <<  std::endl;
  for (unsigned int i=0; i<eta_factor_vector.size(); ++i) {
    EtaCreationFactor *ieta;
    ieta = &eta_factor_vector[i];
    double etaval = ieta->eta_;
    avg_eta += etaval;
    if(etaval>max_eta) {
      max_eta = etaval;
    }
    if(etaval<min_eta) {
      min_eta = etaval;
    }
  }
  avg_eta /= static_cast<double>(eta_factor_vector.size());
  BOOST_LOG_SEV(lg, info) << "birth elements = " << eta_factor_vector.size();
  BOOST_LOG_SEV(lg, info) << "Max eta " << max_eta;
  BOOST_LOG_SEV(lg, info) << "Min eta " << min_eta;
  BOOST_LOG_SEV(lg, info) << "Avg eta " << avg_eta;
}

void CRate::compute_deathfactor()
{
  // death factor
  unsigned int dsum = 0;
// #pragma omp parallel for
  for (unsigned int l = 0; l < gm.vols.nsections; ++l) {
    for (unsigned int q = 0; q < gm.chrgs.nsections; ++q) {
      // WARNING check boundary conditions
      for (unsigned int m = 0; m < gm.vols.nsections-1; ++m) {
        for (unsigned int p = 0; p < gm.chrgs.nsections; ++p) {
          if(rcoag[l][q][m][p] > EPSILONDEATH) {
            dsum++;
            DeathFactor death_factor_aux;
            fill_death(death_factor_aux, dsum, l, q, m, p, rcoag[l][q][m][p]);
            death_factor_vector.push_back(death_factor_aux);
          }
        }
      }
    }
  }

  BOOST_LOG_SEV(lg, info) << "death elements = " << dsum;

  double max_df = death_factor_vector[0].death_;
  double min_df = death_factor_vector[0].death_;
  double avg_df = 0.0;
  // iterate in eta_factor list
//   #pragma omp parallel for shared(max_df, min_df)// schedule(auto)
  for (unsigned int i=0; i<death_factor_vector.size(); ++i) {
    DeathFactor *dfo;
    dfo = &death_factor_vector[i];
    double df = dfo->death_;
    avg_df += df;
    if(df>max_df) {
      max_df = df;
    }
    if(df<min_df) {
      min_df = df;
    }
  }
  avg_df /= static_cast<double>(death_factor_vector.size());
  BOOST_LOG_SEV(lg, info) << "Max death " << max_df;
  BOOST_LOG_SEV(lg, info) << "Min death " << min_df;
  BOOST_LOG_SEV(lg, info) << "Avg death " << avg_df;
}

// Write creation eta to file
void CRate::write_etafactor() {
  
  long int etasize = static_cast<long int>(eta_factor_vector.size());
  bgrid2d etagrid = {{etasize, neta}};
 
  eta2d.resize(etagrid);
  
  unsigned int row=0;
  for(auto ecf : eta_factor_vector) {
      eta2d[row][0] = ecf.id_;
      eta2d[row][1] = ecf.l_;
      eta2d[row][2] = ecf.q_;
      eta2d[row][3] = ecf.m_;
      eta2d[row][4] = ecf.p_;
      eta2d[row][5] = ecf.n_;
      eta2d[row][6] = ecf.r_;
      eta2d[row][7] = ecf.eta_;
      ++row;
  }
  
  int err;
  
  std::string gname = "Electrostatic_interaction";
  std::string dsname = "eta_factor_vector";
  
  err = write_dataset2d_hdf5(h5obj,
                             gname,
                             dsname,
                             eta2d,
                             etasize,
                             neta);

  err = create_dsattrib_hdf5<long int>(h5obj,
                                       gname,
                                       dsname,
                                       "size",
                                       etasize);
}

// Write creation eta to file
void CRate::read_etafactor() {
  
//   long int etasize = static_cast<long int>(eta_factor_vector.size());
  
  int err;
  long int etasize = 0;
  std::string gname = "Electrostatic_interaction";
  std::string dsname = "eta_factor_vector";
 
  err = read_dsattrib_hdf5<long int>(h5obj,
                                     gname,
                                     dsname,
                                     "size",
                                     etasize);

  bgrid2d etagrid = {{etasize, neta}};
  eta2d.resize(etagrid);
  
  err = read_dataset2d_hdf5(h5obj,
                            gname,
                            dsname,
                            eta2d,
                            etasize,
                            neta);
  

  // flatten
  eta_factor_vector.clear();  
  for (unsigned int i = 0; i < etasize; ++i) {
    EtaCreationFactor eta_factor_aux;
    fill_eta(eta_factor_aux,
             eta2d[i][0],
             eta2d[i][1],
             eta2d[i][2],
             eta2d[i][3],
             eta2d[i][4],
             eta2d[i][5],
             eta2d[i][6],
             eta2d[i][7]);
    eta_factor_vector.push_back(eta_factor_aux);
  }

//     for(unsigned int i = 0; i<etasize; ++i) {
//        std::cerr << std::endl;
//       for(unsigned int j = 0; j<neta; ++j) {
// //          array2D[i][j] = varray[i+nrow*j];
//         std::cerr << eta2d[i][j] << '\t';
//       }
//     }
  BOOST_LOG_SEV(lg, info)<< "eta read: " << etasize;
}

// Write death factor to file
void CRate::write_deathfactor() {
 
  long int deathsize = static_cast<long int>(death_factor_vector.size());
  bgrid2d deathgrid = {{deathsize, ndeath}};
//   boost_array2d death2d(deathgrid);
  
  death2d.resize(deathgrid);
  
  unsigned int row=0;
  for(auto ecf : death_factor_vector) {
      death2d[row][0] = ecf.id_;
      death2d[row][1] = ecf.l_;
      death2d[row][2] = ecf.q_;
      death2d[row][3] = ecf.m_;
      death2d[row][4] = ecf.p_;
      death2d[row][5] = ecf.death_;
      ++row;
  }
  
  int err;
  
  std::string gname = "Electrostatic_interaction";
  std::string dsname = "death_factor_vector";
  
  err = write_dataset2d_hdf5(h5obj,
                             gname,
                             dsname,
                             death2d,
                             deathsize,
                             ndeath);

  err = create_dsattrib_hdf5<long int>(h5obj,
                                       gname,
                                       dsname,
                                       "size",
                                       deathsize);
}


// Write death eta to file
void CRate::read_deathfactor() {
  
//   long int etasize = static_cast<long int>(eta_factor_vector.size());
  
  int err;
  long int deathsize = 0;
  
  std::string gname = "Electrostatic_interaction";
  std::string dsname = "death_factor_vector";
 
  err = read_dsattrib_hdf5<long int>(h5obj,
                                     gname,
                                     dsname,
                                     "size",
                                     deathsize);

  bgrid2d deathgrid = {{deathsize, ndeath}};
  death2d.resize(deathgrid);
  
  err = read_dataset2d_hdf5(h5obj,
                            gname,
                            dsname,
                            death2d,
                            deathsize,
                            ndeath);
  
  // flatten
  death_factor_vector.clear();  
  for (unsigned int i = 0; i < deathsize; ++i) {
    DeathFactor death_factor_aux;
    fill_death(death_factor_aux,
               death2d[i][0],
               death2d[i][1],
               death2d[i][2],
               death2d[i][3],
               death2d[i][4],
               death2d[i][5]);
    death_factor_vector.push_back(death_factor_aux);
  }
  
//     for(unsigned int i = 0; i<deathsize; ++i) {
// //        std::cerr << std::endl;
//       for(unsigned int j = 0; j<ndeath; ++j) {
// //          array2D[i][j] = varray[i+nrow*j];
// //         std::cerr << death2d[i][j] << '\t';
//       }
//     }
  BOOST_LOG_SEV(lg, info)<< "death read: " << deathsize;

}


void CRate::write_4d(std::string gname, boost_array4d array4d) {
  
  // Create group
  int err = create_group_hdf5(h5obj, gname);
  int id = 0;
  for (unsigned int l=0; l<gm.vols.nsections; ++l) {
    // iterate in charges particle 1
    for (unsigned int q=0; q<gm.chrgs.nsections; ++q) {
      // iterate in radii particle 2 
      boost_array2d slab(grid);
      for (unsigned int m=0; m<gm.vols.nsections; ++m) {
        // iterate in charges particle 2
        for (unsigned int p=0; p<gm.chrgs.nsections; ++p) {
          slab[m][p] = array4d[l][q][m][p];
        }
      }
      std::string sname = std::to_string(id);
      err = write_dataset2d_hdf5(h5obj,
                                 gname,
                                 sname,
                                 slab,
                                 gm.vols.nsections,
                                 gm.chrgs.nsections);
      err = create_dsattrib_hdf5<unsigned int>(h5obj,
                               gname,
                               sname,
                               "l",
                               l);
      err = create_dsattrib_hdf5<double>(h5obj,
                               gname,
                               sname,
                               "diameter",
                               gm.vols.diameters[l]);
      err = create_dsattrib_hdf5<unsigned int>(h5obj,
                               gname,
                               sname,
                               "q",
                               q);
      err = create_dsattrib_hdf5<double>(h5obj,
                               gname,
                               sname,
                               "charge",
                               gm.chrgs.charges[q]);
      ++id;
    }
  }
}

void CRate::write_efactor() {
  std::string gname = "Enhancement_factor";
  write_4d(gname, efactor);
  int err = create_attrib_hdf5(h5obj, gname,
                               "electhermratio",
                               electhermratio);
}

void CRate::write_efactor_serial() {
  std::string gname = "Enhancement_factor_serial";
  std::string dsname = "Indices";
  //  create_group_hdf5(h5obj, gname);

  write_dset2d_hdf5<boost_short_array2d, short>(h5obj, gname, dsname,
						efindices,
						efindices.shape()[0],
						efindices.shape()[1]);

  std::string edsname = "efactor";
  write_dset_hdf5<darray>(h5obj, gname, edsname,
			  daefactor,
			  efindices.shape()[0]);
   
}

void CRate::write_potentials_serial() {
  std::string cgname = "Contact_potential_serial";
  std::string dsname = "Indices";

  write_dset2d_hdf5<boost_short_array2d, short>(h5obj, cgname, dsname,
						cpindices,
						cpindices.shape()[0],
						cpindices.shape()[1]);
  
  std::string csname = "contact_potential";
  write_dset_hdf5<darray>(h5obj, cgname, csname,
			  dacpotentials,
			  cpindices.shape()[0]);
  

  std::string gname = "Barrier_potential_serial";
  write_dset2d_hdf5<boost_short_array2d, short>(h5obj, gname, dsname,
						bpindices,
						bpindices.shape()[0],
						bpindices.shape()[1]);

  
  std::string bcsname = "contact_potential";
  write_dset_hdf5<darray>(h5obj, gname, bcsname,
			  dabcpotentials,
			  bpindices.shape()[0]);

  std::string bsname = "barrier_potential";
  write_dset_hdf5<darray>(h5obj, gname, bsname,
			  dabpotentials,
			  bpindices.shape()[0]);

  std::string rbsname = "rbarrier";
  write_dset_hdf5<darray>(h5obj, gname, rbsname,
			  darbarriers,
			  bpindices.shape()[0]);

  
}


void CRate::write_potentials() {
  std::string cname = "Contact_potential";
  write_4d(cname, cpotentials);

  std::string bname = "Barrier_potential";
  write_4d(bname, bpotentials);

  std::string rname = "Barrier_location";
  write_4d(rname, rbarriers);  
}

void CRate::write_rcoagulation() {
  std::string gname = "Coagulation_rate";
  write_4d(gname, rcoag);
  int err = create_attrib_hdf5(h5obj, gname,
                               "beta0", beta0);
  
}
