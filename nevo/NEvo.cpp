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
  delete(sol);
  
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

  check = 2;
  // initial conditions, in boost vector_type x( 2 , 1.0 );
  N_Vector xini = NULL;
  xini = N_VNew_Serial(cr.gm.chrgs.nsections);
  NV_Ith_S(xini, 0) = RCONST(2.0);
  NV_Ith_S(xini, 1) = RCONST(1.0);
 
    std::cerr << "\n Initial condition : "
              << NV_Ith_S(xini,0) << '\t' << NV_Ith_S(xini,1) << '\n' << NV_LENGTH_S(xini) << "\n"; 
              
  // Instantiate the solver
  sol = new Solver(this, xini, 0.0);
//  Solver sol(xini, 0.0);
  check = 3;
  sol->print();
  // final time
  realtype tf = RCONST(50.0);
  // delta time
  realtype dt = RCONST(0.01);
//   sol->compute(tf, dt);
//   sol.compute(tf, dt);


  // or solve one step
//   sol.compute_step(0.0, dt);
  
  // free memory for xini
  N_VDestroy_Serial(xini);
  
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
    
    evolve_one_step(ctime);
    
    if(npcount>static_cast<unsigned int>((0.01*nm.tm.tstop)/nm.tm.ndeltat)) {
      err=write_partial_results(ctime);
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
  idens.resize(cr.grid);
  idens[0][cr.gm.chrgs.maxnegative] = nm.cs.indens;

  // auxiliar densities
  pdens.resize(cr.grid);
  pdens = idens;
  ndens.resize(cr.grid);
  
  moments_file = new std::fstream(dirname + nano_filename + "-moments.dat", std::fstream::out);
  *moments_file << "#Time\tNumber\tVolume\tCharge\tmu11\tmu20\tmu02\tmu21\tmu12\tmu30\tmu03";

  moments.resize(10);
  compute_moments(idens);

  // surface growth  
  surface_rate.resize(cr.gm.vols.nsections);
  surface_rate = nm.rs.sgrowth_rate * cr.gm.vols.radii.apply(area_from_radius<double>);

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
}

int NEvo::compute_collisionfreq() {
  double max_efreq = 0.0;
  double max_ifreq = 0.0;
  double min_efreq = 1.0e15;
  double min_ifreq = 1.0e15;

  // From 1.Allen, J. E. Probe theory - the orbital motion approach. Phys. Scr. 45, 497 (1992).
  
  // FIXME optimization use a vector of areas instead of M_PI radii2[l]
  double kte = (2.0/3.0)*pm.es.emean*eCharge;
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
          if ((tfreq[l][q] > efreq[l][q])/*&&(tfreq[l][q] > max_ifreq)*/) {
//              efreq[l][q] = 0.0;
//              tfreq[l][q] = 0.0;
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
//   if(moments[0]!=moments[0]) std::terminate;
  if(!std::isfinite(moments[0])){
    std::cerr << "\n[ee] Particle number is not finite. Terminate.\n";
    BOOST_LOG_SEV(lg, fatal) << "\nParticle number is not finite. Terminate.\n";
//     std::terminate();
  }
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

int NEvo::write_partial_results(double ctime) {

  int err;
  err = write_dataset2d_hdf5(h5obj_nano, "Density", "density"+std::to_string(ctime), ndens,
                             cr.gm.vols.nsections, cr.gm.chrgs.nsections);

  // Current time of simulation
  err = create_attrib_hdf5(h5obj_nano, "Time", "currtime-"+std::to_string(ctime), ctime);

  compute_moments(ndens);
  write_moments(ctime);
}

int NEvo::evolve_one_step(double ctime) {
  int err;
  
  state_type nqdens(cr.gm.chrgs.nsections);
  
//   std::vector<double> nqdens(cr.gm.chrgs.nsections);
  
  double sdt=0.0;
  if(nm.rs.wch == 1) {
//     for(;;) {
//       if(sdt>nm.tm.ndeltat) break;
//       err = compute_explicit_charging(nm.tm.qdeltat);
//   //     std::cerr << sdt;
//       sdt += nm.tm.qdeltat;
//       pdens = ndens;
//     }
    
    
// #pragma omp parallel for
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
  }
  err = advance_nocharging(ctime);
  return err;
}


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

  // BIMODAL
  jnucleation[0][cr.gm.chrgs.maxnegative] = nm.rs.nucleation_rate*1E-3;
  // Warthesen after 10ms no more nucleation, DANGER time splitting and implic
  if(ctime<1e-2) {
//  if(tarea<1.8e-2) {
    wsg = 0.0;
    wnu = 1.0;//1.0
    jnucleation[0][cr.gm.chrgs.maxnegative] = nm.rs.nucleation_rate;
  }

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
  
  if(nm.rs.wco == 1) {
    // explicit
    err = compute_coagulation();

    for (unsigned int l = 0; l < cr.gm.vols.nsections; ++l) {
      for (unsigned int q = 0; q < cr.gm.chrgs.nsections; ++q) {
        ndens[l][q] = pdens[l][q]
                    + wco * nm.tm.ndeltat
                          * kcoagulation[l][q];
        //DANGER WARNING
        if (ndens[l][q]<0.0/*EPSILONETA*/){//EPSILON) {
  //         std::cout << "\n[ee] Negative nanoparticle ndensity (coagulation)";
          std::cerr << "\n[ee] Negative nanoparticle ndensity (coagulation)";
          std::terminate();
          ndens[l][q] = 0.0;
        }
        kcoagulation[l][q] = 0.0;
      }
    }
    
    // update density
    pdens = ndens;
  }
  
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

int NEvo::compute_coagulation() {

#pragma omp parallel
{
#pragma omp for nowait schedule(auto) firstprivate(pdens)// Good
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
#pragma omp for nowait schedule(auto) firstprivate(pdens)// Good
  for (unsigned int j=0; j<cr.death_factor_vector.size(); ++j) {
    //
    DeathFactor *idth;
    idth = &cr.death_factor_vector[j];
//     #pragma omp atomic
    death_vector[idth->l_][idth->q_] += idth->death_
                                      * pdens[idth->m_][idth->p_]
                                      * pdens[idth->l_][idth->q_];
  }

}
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
}

