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

#ifndef CRATE_H
#define CRATE_H

#include <iostream>
#include <string>
#include <functional>
#include <chrono>
#include <ctime>
#include <memory>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

// hdf5 c++ bindings
#include <H5Cpp.h>

#include "utils.h"
#include "log.h"
#include "h5plasma.h"
#include "GridModel.h"
#include "eint.h"
#include "enhancement.h"
#include "eta.h"


/**
  * class CRate
  *
  * Represents the data for coagulation rate calculations
  *
  */
class CRate {
public:
// Constructors/Destructors WARNING, TODO rule of 5
//
  //! Constructor for CRate
  /*!
    @param  dirname_ Directory for output files.
    @param  prefix_filename_ Prefix for output files.
    @param  lg_ Logger instance.
  */
  CRate(std::string dirname_,
        std::string prefix_filename_,
        src::severity_logger< severity_level > lg_)
        : dirname(dirname_),
        prefix_filename(prefix_filename_),
        lg(lg_) {
//     BOOST_LOG_SEV(lg, debug) << "CRate instantiation";
    electhermratio = 0.0;
  }

  //! Constructor for CRate
  /*!
    @param  dirname_ Directory for output files.
    @param  prefix_filename_ Prefix for output files.
    @param  lg_ Logger instance.
  */
  CRate(H5::H5File h5obj_,
        src::severity_logger< severity_level > lg_)
        : h5obj(h5obj_),
        lg(lg_) {
//     BOOST_LOG_SEV(lg, debug) << "CRate instantiation";
  }

  //! Constructor for CRate
  /*!
    @param  dirname_ Directory for output files.
    @param  prefix_filename_ Prefix for output files.
    @param  lg_ Logger instance.
  */
  CRate() {
    
  }
  
// Public methods
//
  //! Open h5 datafile
  /*! 
   * Open the file prefix_filename.h5 in dirname
  */
  int open();

  //! Close h5 datafile
  /*! 
   * Close the file prefix_filename.h5 in dirname
  */
  int close();

  //! Read h5 datafile
  /*! 
   * Read the file prefix_filename.h5 in dirname
  */
  int read();

  //! Read h5 crat results (grid, eta and deathfactor)
  /*! 
   * Read the crat's output file in dirname
  */
  int read_results();

  //! Read and merge h5 file list
  /*! 
   * Read and merge the crat list
  */
  int compute_list(std::vector<std::string>);

  //! Start calculations
  /*! 
   * Start to compute the calculations
  */
  int compute();

  //! Start calculations, need to read pairs first
  /*! 
   * Start to compute the with pairs from a file
  */
  int compute_frompairs();
  
  //! Start calculations WARNING TESTING
  /*! 
   * Testing symmetrization
  */
  int compute_sym();
  
  //! write h5 datafile
  /*! 
   * Write on file prefix_filename.h5 in dirname
  */
  int write();

  //! write h5 datafile from pairs, serialized
  /*! 
   * Write on file prefix_filename.h5 in dirname
  */
  int write_frompairs();
  
  //! write particle pairs to datafile
  /*! 
   * Write 
  */
  int write_pairs();

  //! read particle pairs from datafile
  /*! 
   * Read
  */
  int read_pairs();
  
//
// Public attributes
//
  std::string dirname;                //!< string for directory for output.
  std::string prefix_filename;        //!< string for prefix for output.
  src::severity_logger< severity_level > lg; //!< Logger instance.
  H5::H5File h5obj;                   //!< HDF5 writable object.

  GridModel gm;                       //!< Data model for a grid.

  bgrid2d grid;                       //!< Grid vsecs x qsecs
    
  bgrid4d grid4;                      //!< Grid vsecs x qsecs x vsecs x qsecs

  //!< Container for EtaCreationFactor.
  std::vector<EtaCreationFactor> eta_factor_vector;

  //!< Container for DeathFactor.
  std::vector<DeathFactor> death_factor_vector;
//
private:
// Private methods
//
  //! Bind method to efactor function
  /*!
   *  Method to efactor funct
  */
  void bind_efactorfunc();

  //! Compute efactor
  /*! 
    Compute enhancement factor using efactorfunc
    @param  r1 radius of particle 1
    @param  r2 radius of particle 2
    @param  q1 charge of particle 1
    @param  q2 charge of particle 2
  */
  double compute_efactor(double r1, double r2,
                     double q1, double q2) {
    return (this->*efactorfunc)(r1, r2,
                                q1, q2);

  }

  //! Compute efactor using mpc method
  /*!
    Compute enhancement factor using multipolar coefficient expansion
    @param  r1 radius of particle 1
    @param  r2 radius of particle 2
    @param  q1 charge of particle 1
    @param  q2 charge of particle 2
  */
  double compute_efactor_mpc(double r1, double r2,
                             double q1, double q2) {

    /* return eint::efactor_mpc(r1, r2, */
    /*                          q1, q2, */
    /*                          gm.einter.dconstant, */
    /*                          gm.gsys.temperature, */
    /*                          gm.einter.terms, */
    /*                          gm.einter.terms); */
    return 1.0;
  }

  //! Compute efactor using ipa method
  /*!
    Compute enhancement factor using image potential approximation expansion
    @param  r1 radius of particle 1
    @param  r2 radius of particle 2
    @param  q1 charge of particle 1
    @param  q2 charge of particle 2
  */
  double compute_efactor_ipa(double r1, double r2,
                             double q1, double q2) {

    /* return eint::efactor_ipa(r1, r2, */
    /*                          q1, q2, */
    /*                          gm.einter.dconstant, */
    /*                          gm.gsys.temperature, */
    /* 			     gm.einter.terms, */
    /*                          gm.einter.terms); */
    return 1.0;
  }

  //! Compute efactor using Coulomb method
  /*! 
    Compute enhancement factor using Coulomb potential
    @param  r1 radius of particle 1
    @param  r2 radius of particle 2
    @param  q1 charge of particle 1
    @param  q2 charge of particle 2
  */
  double compute_efactor_coul(double r1, double r2,
                              double q1, double q2) {
    if(q1*q2>0) {
       return exp(-electhermratio*q1*q2/(r1+r2));
    }
    else {
       return 1.0-electhermratio*q1*q2/(r1+r2);
    }
  }

  //! Compute efactor in the grid vols x chrgs
  /*! 
    Compute enhancement in the grid vols x chrgs
  */
  void compute_efactor_grid();

  //! Compute efactor in the grid vols x chrgs
  /*! 
    Compute enhancement in the grid vols x chrgs symmetric
  */
  void compute_efactor_grid_sym();

  double beta_free(const double vol1, const double vol2,
                   const double beta0);

  //! Compute coagulation rate
  /*! 
    Compute coagulation rate
  */
  void compute_rcoagulation();

  //! Compute etafactor
  /*! 
    Compute etafactor
  */
  void compute_etafactor();

  //! Compute etafactor volume and charge pivoting
  /*! 
    Compute etafactor two components volume and charge
  */
  void compute_etafactor2C();
  
  //! Write etafactor
  /*! 
    Write etafactor
  */
  void write_etafactor();

  //! Read etafactor
  /*! 
    Read etafactor
  */
  void read_etafactor();
  
  //! Compute deathfactor
  /*! 
    Compute deathfactor
  */
  void compute_deathfactor();


  //! Compute deathfactor volume and charge pivoting
  /*! 
    Compute deathfactor two components volume and charge
  */
  void compute_deathfactor2C();
  
  //! Write deathfactor
  /*! 
    Write deathfactor
  */
  void write_deathfactor();

  //! Read deathfactor
  /*! 
    Read deathfactor
  */
  void read_deathfactor();

  //! Write boost_array4d
  /*! 
    Write boost_array4d
  */
  void write_4d(std::string gname, boost_array4d array4d);

  //! Write enhancement factor
  /*! 
    Write enhancement factor
  */
  void write_efactor();

  //! Write enhancement factor
  /*! 
    Write enhancement factor
  */
  void write_efactor_serial();

  //! Write potentials
  /*! 
    Writes contact and barrier potentials
  */
  void write_potentials();

  //! Write potentials
  /*! 
    Writes contact and barrier potentials
  */
  void write_potentials_serial();
  
  //! Write coagulation rate
  /*! 
    Write coagulation rate
  */
  void write_rcoagulation();

//
// Private attributes
// 

  boost_array4d efactor;              //!< Array for enhancement factor.
  boost_array4d cpotentials;          //!< Array for contact potentials.
  boost_array4d bpotentials;          //!< Array for potential barriers.
  boost_array4d rbarriers;            //!< Array for potential barriers.

  boost_short_array2d efindices;    //!< Array for enhancement factor indices.
  boost_short_array2d cpindices;    //!< Array for contact potential indices.
  boost_short_array2d bpindices;    //!< Array for rbarrier and barrier potential indices.
    
  darray daefactor;
  darray dacpotentials;
  darray dabpotentials;
  darray dabcpotentials;
  darray darbarriers;
  
  /* boost_array4d_ref efactor_ref;     //!< Reference to efactor */
  boost_array4d rcoag;                //!< Coagulation rate.

  double electhermratio;          //!< Electrostatic energy to thermal ratio.

  //!< Type definition to pointer to function.
  typedef double ( CRate::*efactorfunc_t ) ( double, double, double, double );  

  efactorfunc_t efactorfunc;      //!< Pointer to function.

  double beta0;                   //!< Coagulation rate beta0 prefactor.

  //!< Container for EtaFactor.
  boost_array2d eta2d;

  //!< Container for DeathFactor.
  boost_array2d death2d;

  static const unsigned int neta = 8; //!< number of eta indices.
  static const unsigned int ndeath = 6;//!< number of death indices.

  std::vector<std::string> sfilelist;
};

#endif // CRATE_H
