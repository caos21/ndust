/*
 * Copyright 2016 <Benjamin Santos> <caos21@gmail.com>
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

#include "GridModel.h"

int GridModel::read() {

  int err;

  BOOST_LOG_SEV(lg, info) << "Reading file. Id: "
                          << h5obj.getId();

  // description attributes
  err = read_description();

  // gridsystem attributes
  err = read_gridsystem();

  // vsections attributes
  err = read_vsections();

  // qsections attributes
  err = read_qsections();

  // Einteraction attributes
  err = read_einteraction();
//   if(err!=0) {
//     BOOST_LOG_SEV(lg, error) << "h5 I/O error closiing file. Terminate.";
//     std::terminate();
//   }
//   BOOST_LOG_SEV(lg, info) << "Closed file Id: " << id;
  
  return 0;
}

//
// ----------------------------- Description ----------------------------------
int GridModel::read_description() {
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
// ------------------------------- GridModel ----------------------------------
int GridModel::read_gridsystem() {
  int err = 0;
  BOOST_LOG_SEV(lg, info) << "Grid system: ";

  err = read_attrib_hdf5<double>(h5obj, "/Grid_system", "Temperature",
                         gsys.temperature);
  BOOST_LOG_SEV(lg, info) << "--> Temperature: " << gsys.temperature;

  err = read_attrib_hdf5<double>(h5obj, "/Grid_system", "Mass_density",
                         gsys.nmdensity);
  BOOST_LOG_SEV(lg, info) << "--> Nanoparticle mass density: "
                          << gsys.nmdensity;

  return 0;
}

//
// ------------------------------- VSections ----------------------------------
int GridModel::read_vsections() {
  int err = 0;
  BOOST_LOG_SEV(lg, info) << "Volume sections: ";

  err = read_attrib_hdf5<unsigned int>(h5obj, "/Volume_sections", "NSections",
                         vols.nsections);
  BOOST_LOG_SEV(lg, info) << "--> Number of sections: " << vols.nsections;

  err = read_attrib_hdf5<double>(h5obj, "/Volume_sections", "Base",
                         vols.base);
  BOOST_LOG_SEV(lg, info) << "--> Base: " << vols.base;

  err = read_attrib_hdf5<double>(h5obj, "/Volume_sections", "Power",
                         vols.power);
  BOOST_LOG_SEV(lg, info) << "--> Power: " << vols.power;

  err = read_attrib_hdf5<double>(h5obj, "/Volume_sections", "Min_interface",
                         vols.miniface);  
  BOOST_LOG_SEV(lg, info) << "--> Minimum volume at interface: "
                          << vols.miniface;
  
  err = read_dataset_hdf5<darray>(h5obj, vols.interfaces, "/Volume_sections", "Interfaces");
  BOOST_LOG_SEV(lg, info) << "--> Interfaces: Read " << vols.interfaces.size() << " elements";
  
  err = read_dataset_hdf5<darray>(h5obj, vols.volumes, "/Volume_sections", "Volumes");
  BOOST_LOG_SEV(lg, info) << "--> Volumes: Read " << vols.volumes.size() << " elements";
  
  err = read_dataset_hdf5<darray>(h5obj, vols.radii, "/Volume_sections", "Radii");
  BOOST_LOG_SEV(lg, info) << "--> Radii: Read " << vols.radii.size() << " elements";
  
  err = read_dataset_hdf5<darray>(h5obj, vols.diameters, "/Volume_sections", "Diameters");
  BOOST_LOG_SEV(lg, info) << "--> Diameters: Read " << vols.diameters.size() << " elements";

  BOOST_LOG_SEV(lg, info) << "--> Minimum diameter: " << vols.diameters[0];
  BOOST_LOG_SEV(lg, info) << "--> Maximum diameter: " << vols.diameters[vols.nsections-1];
  
  return 0;
}

//
// ------------------------------- QSections ----------------------------------
int GridModel::read_qsections() {
  int err = 0;
  BOOST_LOG_SEV(lg, info) << "Charge sections: ";
  
  err = read_attrib_hdf5<unsigned int>(h5obj, "/Charge_sections", "NSections",
                         chrgs.nsections);
  BOOST_LOG_SEV(lg, info) << "--> Number of sections: " << chrgs.nsections;

  err = read_attrib_hdf5<unsigned int>(h5obj, "/Charge_sections", "Max_negative",
                         chrgs.maxnegative);  
  BOOST_LOG_SEV(lg, info) << "--> Maximum negative charge: "
                          << chrgs.maxnegative;

  err = read_attrib_hdf5<unsigned int>(h5obj, "/Charge_sections", "Max_positive",
                         chrgs.maxpositive);  
  BOOST_LOG_SEV(lg, info) << "--> Maximum positive charge: "
                          << chrgs.maxpositive;

  err = read_dataset_hdf5<darray>(h5obj, chrgs.charges, "/Charge_sections", "Charges");
  BOOST_LOG_SEV(lg, info) << "--> Charges: Read " << chrgs.charges.size() << " elements";
  
  return 0;
}

//
// ----------------------------- Einteraction ---------------------------------
int GridModel::read_einteraction() {
  int err = 0;
  BOOST_LOG_SEV(lg, info) << "Electrostatic interaction: ";

  err = read_attrib_hdf5<double>(h5obj, "/Electrostatic_interaction", "Multiplier",
                         einter.multiplier);
  BOOST_LOG_SEV(lg, info) << "--> Multiplier: "
                          << einter.multiplier;
                          
  err = read_attrib_hdf5<double>(h5obj, "/Electrostatic_interaction", "Dielectric_constant",
                         einter.dconstant);  
  BOOST_LOG_SEV(lg, info) << "--> Dielectric constant: "
                          << einter.dconstant;
                          
  err = read_attrib_hdf5<unsigned int>(h5obj, "/Electrostatic_interaction", "Method",
                         einter.method);
  BOOST_LOG_SEV(lg, info) << "--> Method: " << einter.method;

  // terms make sense only for MPC, method 0
  if(einter.method == 0) {
    err = read_attrib_hdf5<unsigned int>(h5obj, "/Electrostatic_interaction", "Terms",
                          einter.terms);
    BOOST_LOG_SEV(lg, info) << "--> Terms for MPC: " << einter.terms;
  }

  return 0;
}
