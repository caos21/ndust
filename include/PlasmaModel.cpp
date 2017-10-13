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

#include "PlasmaModel.h"
  
int PlasmaModel::read() {

  int err;
  
  BOOST_LOG_SEV(lg, info) << "Reading file. Id: "
                          << h5obj.getId();
  
  BOOST_LOG_SEV(lg, info) << "--> Description: " << desc.description;
  // description attributes
  err = read_description();

  // read parameters
  err = read_parameters();
  
  // electrons attributes
  err = read_electrons();

  // ions attributes
  err = read_ions();

  // meta attributes
  err = read_metastables();
  
  return 0;
}

//
// ----------------------------- Description ----------------------------------
int PlasmaModel::read_description() {
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
// ------------------------------- PlasmaModel ----------------------------------
int PlasmaModel::read_parameters() {
  int err = 0;
  BOOST_LOG_SEV(lg, info) << "Plasma parameters: ";
  
  err = read_attrib_hdf5<double>(h5obj, "/Parameters", "length",
                         pars.length);
  
  BOOST_LOG_SEV(lg, info) << "--> Length: " << pars.length;


  err = read_attrib_hdf5<double>(h5obj, "/Parameters",
                                 "neutral_density",
                                 pars.neutral_density);
  
  BOOST_LOG_SEV(lg, info) << "--> Neutral density: " << pars.neutral_density;


  err = read_attrib_hdf5<double>(h5obj, "/Parameters", "pressure",
                                pars.pressure);
  
  BOOST_LOG_SEV(lg, info) << "--> Pressure: " << pars.pressure;  

  err = read_attrib_hdf5<double>(h5obj, "/Parameters", "temperature",
                                pars.temperature);
  
  BOOST_LOG_SEV(lg, info) << "--> Temperature: " << pars.temperature;

  return 0;
}

//
// ------------------------------- electrons ----------------------------------
int PlasmaModel::read_electrons() {
  int err = 0;
  BOOST_LOG_SEV(lg, info) << "Electrons: ";
  
  err = read_attrib_hdf5<double>(h5obj, "/Electrons", "emean",
                                es.emean);
  
  BOOST_LOG_SEV(lg, info) << "--> emean: " << es.emean;


  err = read_attrib_hdf5<double>(h5obj, "/Electrons", "ne",
                                es.ne);
  
  BOOST_LOG_SEV(lg, info) << "--> ne: " << es.ne;

  return 0;
}

//
// ------------------------------- ions ---------------------------------------
int PlasmaModel::read_ions() {
  int err = 0;
  BOOST_LOG_SEV(lg, info) << "Ions: ";
  
  err = read_attrib_hdf5<double>(h5obj, "/Ions", "itemp",
                                 is.itemp);
  
  BOOST_LOG_SEV(lg, info) << "--> ion temperature: " << is.itemp;


  err = read_attrib_hdf5<double>(h5obj, "/Ions", "imass",
                                 is.imass);
  
  BOOST_LOG_SEV(lg, info) << "--> ion mass: " << is.imass;


  err = read_attrib_hdf5<double>(h5obj, "/Ions", "ni",
                                 is.ni);
  
  BOOST_LOG_SEV(lg, info) << "--> ni: " << is.ni;

  return 0;
}

//
// ------------------------------- Metastables --------------------------------
int PlasmaModel::read_metastables() {
  int err = 0;
  BOOST_LOG_SEV(lg, info) << "Metastables: ";

  err = read_attrib_hdf5<double>(h5obj, "/Metastables", "nm",
                                 ms.nm);
  
  BOOST_LOG_SEV(lg, info) << "--> nm: " << ms.nm;

  return 0;
}
