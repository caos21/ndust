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

#include "NanoModel.h"

int NanoModel::read() {

  int err;

  BOOST_LOG_SEV(lg, info) << "Reading file. Id: "
                          << h5obj.getId();

  BOOST_LOG_SEV(lg, info) << "--> Description: " << desc.description;
  // description attributes
  err = read_description();

  // read nanoparticles
  err = read_nanoparticles();

  // read rates attributes
  err = read_rates();

  // read density attributes
  err = read_density();

  // read time attributes
  err = read_time();

  return 0;
}

//
// ----------------------------- Description ----------------------------------
int NanoModel::read_description() {
  int err = 0;
  BOOST_LOG_SEV(lg, info) << "Sysinfo and description: ";

  err = read_str_attrib_hdf5(h5obj, "/", "Nodename", desc.nodename);
  BOOST_LOG_SEV(lg, info) << "--> Nodename: " << desc.nodename;

  err = read_str_attrib_hdf5(h5obj, "/", "Sysname", desc.sysname);
  BOOST_LOG_SEV(lg, info) << "--> Sysname: " << desc.sysname;

  err = read_str_attrib_hdf5(h5obj, "/", "Machine", desc.machine);
  BOOST_LOG_SEV(lg, info) << "--> Machine: " << desc.machine;

  err = read_str_attrib_hdf5(h5obj, "/", "Release", desc.release);
  BOOST_LOG_SEV(lg, info) << "--> Release: " << desc.release;

  err = read_str_attrib_hdf5(h5obj, "/", "Version", desc.version);
  BOOST_LOG_SEV(lg, info) << "--> Version: " << desc.version;

  err = read_str_attrib_hdf5(h5obj, "/", "Timestamp", desc.timestamp);
  BOOST_LOG_SEV(lg, info) << "--> Timestamp: " << desc.timestamp;

  err = read_str_attrib_hdf5(h5obj, "/", "Description", desc.description);
  BOOST_LOG_SEV(lg, info) << "--> Description: " << desc.description;

  return 0;
}

//
// ----------------------------- nanoparticles --------------------------------
int NanoModel::read_nanoparticles() {
  int err = 0;
  BOOST_LOG_SEV(lg, info) << "Nanoparticle parameters: ";

  err = read_attrib_hdf5<double>(h5obj, "/Nanoparticles", "accfactor",
                                 nano.accfactor);

  BOOST_LOG_SEV(lg, info) << "--> Accumulation factor: " << nano.accfactor;

  err = read_attrib_hdf5<double>(h5obj, "/Nanoparticles", "eaffinity",
                                 nano.eaffinity);

  BOOST_LOG_SEV(lg, info) << "--> Electron affinity: " << nano.eaffinity;

  err = read_attrib_hdf5<int>(h5obj, "/Nanoparticles", "tunnel",
                              nano.tunnel);

  BOOST_LOG_SEV(lg, info) << "--> Tunnel effect: "
                          << (nano.tunnel==1 ? "yes" : "no");

  return 0;
}

//
// ---------------------------------- rates -----------------------------------
int NanoModel::read_rates() {
  int err = 0;
  BOOST_LOG_SEV(lg, info) << "Rates: ";

  err = read_attrib_hdf5<double>(h5obj, "/Rates", "nucleation_rate",
                                 rs.nucleation_rate);

  BOOST_LOG_SEV(lg, info) << "--> Nucleation rate: " << rs.nucleation_rate;


  err = read_attrib_hdf5<double>(h5obj, "/Rates", "sgrowth_rate",
                                 rs.sgrowth_rate);

  BOOST_LOG_SEV(lg, info) << "--> Surface growth rate: " << rs.sgrowth_rate;

  err = read_attrib_hdf5<int>(h5obj, "/Rates", "wch",
                              rs.wch);

  BOOST_LOG_SEV(lg, info) << "--> Nanoparticle charging: "
                          << (rs.wch==1 ? "yes" : "no");

  err = read_attrib_hdf5<int>(h5obj, "/Rates", "wco",
                              rs.wco);

  BOOST_LOG_SEV(lg, info) << "--> Coagulation: "
                          << (rs.wco==1 ? "yes" : "no");

  err = read_attrib_hdf5<int>(h5obj, "/Rates", "wnu",
                              rs.wnu);

  BOOST_LOG_SEV(lg, info) << "--> Nucleation: "
                          << (rs.wnu==1 ? "yes" : "no");

  err = read_attrib_hdf5<int>(h5obj, "/Rates", "wsg",
                              rs.wsg);

  BOOST_LOG_SEV(lg, info) << "--> Surface growth: "
                          << (rs.wsg==1 ? "yes" : "no");

  err = read_attrib_hdf5<int>(h5obj, "/Rates", "wsih4",
                              rs.wsih4);

  BOOST_LOG_SEV(lg, info) << "--> SiH4 coupling: "
                          << (rs.wsih4==1 ? "yes" : "no");

  err = read_attrib_hdf5<double>(h5obj, "/Rates", "sih4ratio",
                                 rs.sih4ratio);

  BOOST_LOG_SEV(lg, info) << "--> SiH4 to gas ratio: " << rs.sih4ratio;
  
  err = read_attrib_hdf5<int>(h5obj, "/Rates", "sih4nmol",
                              rs.sih4nmol);

  BOOST_LOG_SEV(lg, info) << "--> SiH4 molecules per nucleated nanoparticle: "
                          << rs.sih4nmol;

  err = read_attrib_hdf5<double>(h5obj, "/Rates", "sih4mass",
				 rs.sih4mass);

  BOOST_LOG_SEV(lg, info) << "--> SiH4 mass: " << rs.sih4mass;
  return 0;
}

//
// -------------------------------- Density -----------------------------------
int NanoModel::read_density() {
  int err = 0;
  std::string sgroup("/Density");

  BOOST_LOG_SEV(lg, info) << "Density parameters: ";

  err = read_attrib_hdf5<double>(h5obj, sgroup, "indens",
                                 ds.indens);

  BOOST_LOG_SEV(lg, info) << "--> Nanoparticle initial density: " << ds.indens;


  err = read_attrib_hdf5<double>(h5obj, sgroup, "qtol",
                                 ds.qtol);

  BOOST_LOG_SEV(lg, info) << "--> Nanoparticle qtol: " << ds.qtol;

  err = read_attrib_hdf5<int>(h5obj, sgroup, "distribution",
                              ds.distribution);

  BOOST_LOG_SEV(lg, info) << "--> Distribution type: "
                          << (ds.distribution == 0 ? "delta" :
                              (ds.distribution == 1 ? "step" : "Gaussian"));

  err = read_attrib_hdf5<unsigned int>(h5obj, sgroup, "peakpos", ds.peakpos);

  BOOST_LOG_SEV(lg, info) << "--> Peak in section #: " << ds.peakpos;

  if(ds.distribution > 0) {
    err = read_attrib_hdf5<unsigned int>(h5obj, sgroup, "width", ds.width);

    BOOST_LOG_SEV(lg, info) << "--> Distribution width: " << ds.width;
  }

  unsigned int wchargewidth = 0;
  err = read_attrib_hdf5<unsigned int>(h5obj, sgroup, "chargewidth", wchargewidth);
  if (wchargewidth>0) {
    ds.chargewidth = true;
    BOOST_LOG_SEV(lg, info) << "--> Distribution charge width: " << ds.chargewidth;
    err = read_attrib_hdf5<int>(h5obj, sgroup, "chargenegwidth", ds.chargenegwidth);
    BOOST_LOG_SEV(lg, info) << "    + Maximum negative charge width: " << ds.chargenegwidth;
    err = read_attrib_hdf5<int>(h5obj, sgroup, "chargeposwidth", ds.chargeposwidth);
    BOOST_LOG_SEV(lg, info) << "    + Maximum positive charge width: " << ds.chargeposwidth;
  }
  else {
    ds.chargewidth = false;
    ds.chargenegwidth = 0;
    ds.chargeposwidth = 0;
    BOOST_LOG_SEV(lg, info) << "--> Distribution charge width: " << ds.chargewidth;
  }  

  return 0;
}

//
// --------------------------------- Time -------------------------------------
int NanoModel::read_time() {
  int err = 0;
  BOOST_LOG_SEV(lg, info) << "Time parameters: ";

  err = read_attrib_hdf5<double>(h5obj, "/Time", "ndeltat",
                                 tm.ndeltat);

  BOOST_LOG_SEV(lg, info) << "--> Nanoparticle delta t: " << tm.ndeltat;


  err = read_attrib_hdf5<double>(h5obj, "/Time", "qdeltat",
                                 tm.qdeltat);

  BOOST_LOG_SEV(lg, info) << "--> Charging delta t: " << tm.qdeltat;

  err = read_attrib_hdf5<double>(h5obj, "/Time", "tstop",
                                 tm.tstop);

  BOOST_LOG_SEV(lg, info) << "--> Stop time: " << tm.tstop;


  return 0;
}
