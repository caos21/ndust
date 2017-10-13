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

#ifndef PLASMAMODEL_H
#define PLASMAMODEL_H

#include <iostream>
#include <string>

// hdf5 c++ bindings
#include <H5Cpp.h>

#include "log.h"
#include "array.h"
#include "h5plasma.h"
#include "description.h"


/**
  * struct parameters
  *
  * Stores some system parameters
  *
  */
struct parameters {
  double length;                    //!< Reactor length
  double neutral_density;           //!< Neutral gas density
  bool pfixed;                      //!< Plasma fixed conditions
  double pressure;                  //!< Reactor pressure
  double temperature;               //!< Reactor temperature
};

/**
  * struct electrons
  *
  * Stores electrons density 
  *
  */
struct electrons {
  double emean;                     //!< Electron mean energy
  double ne;                        //!< Electron density
};

/**
  * struct ions
  *
  * Stores ions density 
  *
  */
struct ions {
  double itemp;                     //!< Ion temperature
  double imass;                     //!< Ion mass
  double ni;                        //!< Ion density  
};

/**
  * struct metastables
  *
  * Stores metastables density 
  *
  */
struct metastables {
  double nm;                        //!< Electron density
};

/**
  * class PlasmaModel
  *
  * Represents the data model for a plasma calculation
  *
  */
class PlasmaModel {
public:
// Constructors/Destructors
//
  //! Constructor for PlasmaModel
  /*!
    @param  h5obj HDF5 writable object.
    @param  lg_ Logger instance.
  */
  PlasmaModel(H5::H5File h5obj_, src::severity_logger< severity_level > lg_)
        : h5obj(h5obj_), lg(lg_)  {
//     BOOST_LOG_SEV(lg, debug) << "PlasmaModel instantiation";
  }

  PlasmaModel() { }

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

  description desc;                   //!< Data and file description.
  parameters pars;                    //!< System parameters.
  electrons es;                       //!< Electrons.
  ions is;                            //!< Ions.
  metastables ms;                //!< Electrostatic interaction.
  //
private:
// Private methods
//
  //! Read h5 datafile description
  /*! 
   * Read the datafile description text and sysinfo attributes
  */
  int read_description();

  //! Read h5 datafile parameters
  /*! 
   * Read the datafile parameters attributes
  */
  int read_parameters();
  
  //! Read h5 datafile electrons
  /*! 
   * Read the datafile electrons attributes
  */
  int read_electrons();

  //! Read h5 datafile ions
  /*! 
   * Read the datafile ions attributes
  */
  int read_ions();

  //! Read h5 datafile metastables
  /*! 
   * Read the datafile metastables attributes
  */
  int read_metastables();

};

#endif
