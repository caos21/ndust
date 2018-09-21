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

#ifndef GRIDMODEL_H
#define GRIDMODEL_H

#include <iostream>
#include <string>

// hdf5 c++ bindings
#include <H5Cpp.h>

#include "log.h"
#include "array.h"
#include "h5plasma.h"
#include "description.h"

/**
  * struct gridsystem
  *
  * Stores some system parameters
  *
  */
struct gridsystem {
  double temperature;               //!< Nanoparticle temperature
  double nmdensity;                 //!< Nanoparticle mass density
};

/**
  * struct vsections
  *
  * Stores volume section information
  *
  */
struct vsections {
  unsigned int nsections;           //!< Number of volume sections
  double miniface;                  //!< Minimum volume at interface
  double base;                      //!< Base
  double power;                     //!< Power
  darray interfaces;                //!< Array of interfaces
  darray volumes;                   //!< Array of volumes
  darray radii;                     //!< Array of radii
  darray diameters;                 //!< Array of diameters
};

/**
  * struct qsections
  *
  * Stores charge section information
  *
  */
struct qsections {
  unsigned int nsections;           //!< Number of charge sections
  unsigned int maxpositive;         //!< Maximum positive charge
  unsigned int maxnegative;         //!< Maximum negative charge
  darray charges;                   //!< Array of charges
};

/**
  * struct einteraction
  *
  * Stores electronic interaction information
  *
  */
struct einteraction {
  double multiplier;                //!< Multiplier for interaction effect
  double dconstant;                 //!< Dielectric constant for nanoparticles
  unsigned int method;              //!< Method used to compute interaction
  unsigned int terms;               //!< Number of terms of MPC expansion
  double hamaker;                   //!< Hamaker constant
  double vdw_radius;                //!< Van der Waals radius or cutoff for potential
};
/**
  * class GridModel
  *
  * Represents the data model for a grid calculation
  *
  */
class GridModel {
public:
// Constructors/Destructors
//
  //! Constructor for GridModel
  /*!
    @param  h5obj HDF5 writable object.
    @param  lg_ Logger instance.
  */
  GridModel(H5::H5File h5obj_, src::severity_logger< severity_level > lg_)
        : h5obj(h5obj_), lg(lg_)  {
//     BOOST_LOG_SEV(lg, debug) << "GridModel instantiation";
  }

  GridModel() { }

// Public methods
//
  //! Read h5 datafile
  /*! 
   * Read the grid model
  */
  int read();
  
// Public attributes
//
  H5::H5File h5obj;                   //!< HDF5 writable object.
  src::severity_logger< severity_level > lg; //!< Logger instance.

  description desc;                   //!< Data and file description.
  gridsystem gsys;                    //!< Grid system attributes.
  vsections vols;                     //!< Volume sections.
  qsections chrgs;                    //!< Charge sections.
  einteraction einter;                //!< Electrostatic interaction.
  //
private:
// Private methods
//
  //! Read h5 datafile description
  /*! 
   * Read the datafile description text and sysinfo attributes
  */
  int read_description();

  //! Read h5 datafile gridsystem
  /*! 
   * Read the datafile gridsystem attributes
  */
  int read_gridsystem();
  
  //! Read h5 datafile volume sections
  /*! 
   * Read the datafile volume sections and attributes
  */
  int read_vsections();

  //! Read h5 datafile charge sections
  /*! 
   * Read the datafile charge sections and attributes
  */
  int read_qsections();
  
  //! Read h5 datafile electrostatic interaction
  /*! 
   * Read the electrostatic interaction attributes
  */
  int read_einteraction();
};

#endif
