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

#ifndef NANOMODEL_H
#define NANOMODEL_H

#include <iostream>
#include <string>

// hdf5 c++ bindings
#include <H5Cpp.h>

#include "log.h"
#include "array.h"
#include "h5plasma.h"
#include "description.h"


/**
  * struct nanoparticles
  *
  * Stores some information about nanoparticles
  *
  */
struct nanoparticles {
  double accfactor;                 //!< Accomodation factor
  double eaffinity;                 //!< Electron affinity
  int tunnel;                       //!< Tunnel effect switch
};

/**
  * struct rates
  *
  * Stores rates
  *
  */
struct rates {
  double nucleation_rate;           //!< Nucleation rate in 
  double sgrowth_rate;              //!< Surface Growth
  int wch;                          //!< Charging switch
  int wco;                          //!< Coagulation switch
  int wnu;                          //!< Nucleation switch
  int wsg;                          //!< Surface growth switch
  int wsih4;                        //!< SiH4 coupled rates switch
  double sih4ratio;                 //!< SiH4 to ratio
  int sih4nmol;                     //!< Number of SiH4 per nucleated particle
  double sih4mass;                  //!< SiH4 mass
};

/**
  * struct time
  *
  * Stores time parameters 
  *
  */
struct times {
  double ndeltat;                   //!< Nanoparticle delta t
  double qdeltat;                   //!< Charging delta t
  double tstop;                     //!< Time of simulation
};

/**
  * struct density
  *
  * Stores density parameters 
  *
  */
struct density {
  double indens;                   //!< Initial nanoparticle density
  double qtol;                     //!< Charging qtol
  int distribution;                /*!< Type of distribution:
                                        0->delta, 1->step, 2->Gaussian*/
  unsigned int peakpos;            /*!< Number section for peak of initial
                                        Nanoparticle distribution*/
  unsigned int width;              //!< Number section for peak of initial

  bool chargewidth;               //!< Charge width for distribution
  int chargenegwidth;             //!< Max negative charge
  int chargeposwidth;             //!< Max positive charge
  
};


/**
  * class NanoModel
  *
  * Represents the data model for a Nanoparticle calculation
  *
  */
class NanoModel {
public:
// Constructors/Destructors
//
  //! Constructor for NanoModel
  /*!
   * @param  h5obj HDF5 writable object.
   * @param  lg_ Logger instance.
  */
  NanoModel(H5::H5File h5obj_, src::severity_logger< severity_level > lg_)
        : h5obj(h5obj_), lg(lg_)  {
//     BOOST_LOG_SEV(lg, debug) << "NanoModel instantiation";
  }

  NanoModel() { }

// Public methods
//
  //! Read h5 datafile
  /*!
   * Read the plasma model
  */
  int read();

// Public attributes
//
  H5::H5File h5obj;                   //!< HDF5 writable object.
  src::severity_logger< severity_level > lg; //!< Logger instance.

  description desc;               //!< Data and file description.
  nanoparticles nano;             //!< Nanoparticle parameters.
  rates rs;                       //!< Growth rate information.
  density ds;                     //!< Initial density parameters.
  times tm;                       //!< Time parameters.
  //
private:
// Private methods
//
  //! Read h5 datafile description
  /*!
   * Read the datafile description text and sysinfo attributes
  */
  int read_description();

  //! Read h5 datafile nanoparticle parameters
  /*!
   * Read the datafile nanoparticle attributes
  */
  int read_nanoparticles();

  //! Read h5 datafile rates
  /*!
   * Read the datafile rates attributes
  */
  int read_rates();

  //! Read h5 datafile density group
  /*!
   * Read the datafile density attributes
  */
  int read_density();

  //! Read h5 datafile time
  /*!
   * Read the datafile time attributes
  */
  int read_time();

};

#endif
