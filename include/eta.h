/*
 * Copyright 2016 Benjamin Santos <caos21@gmail.com>
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

#ifndef ETA_H
#define ETA_H

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <H5Cpp.h>

#include <sstream>
#include <fstream>
//! Eta factor
/*! Represents the eta factor
 *
 */
class EtaCreationFactor {
public:

  EtaCreationFactor() {
    // call fill to set 0 all the fields
    fill_eta(*this, 0, 0, 0, 0, 0, 0, 0, 0.0);
  }
  unsigned int id_;           //!< Identifier
  short l_;            //!< Index l of current volume section
  short q_;            //!< Index q of current charge section
  short m_;            //!< Index m of coalescing volume
  short p_;            //!< Index p of coalescing charge
  short n_;            //!< Index n of coalescing volume
  short r_;            //!< Index r of coalescing volume
  double eta_;                //!< Value for \f$\eta_{mp,nr}(v_l, Q_q)\f$

  friend void fill_eta(EtaCreationFactor &ecf,
                       const unsigned int id,
                       const short l,
                       const short q,
                       const short m,
                       const short p,
                       const short n,
                       const short r,
                       const double eta);

  friend void etafactor_tostream(const EtaCreationFactor &ecf,
                                 std::fstream& etafactor_file);

  friend void etafactor_fromstream(EtaCreationFactor &ecf,
                                   std::istringstream &ss);
};

inline void fill_eta(EtaCreationFactor &ecf,
                     const unsigned int id,
                     const short l,
                     const short q,
                     const short m,
                     const short p,
                     const short n,
                     const short r,
                     const double eta) {
  ecf.id_=id;
  ecf.l_=l;
  ecf.q_=q;
  ecf.m_=m;
  ecf.p_=p;
  ecf.n_=n;
  ecf.r_=r;
  ecf.eta_=eta;
}

inline
void etafactor_tostream(const EtaCreationFactor &ecf, std::fstream& etafactor_file) {
  etafactor_file << ecf.id_ << '\t'
                 << ecf.l_ << '\t'
                 << ecf.q_ << '\t'
                 << ecf.m_ << '\t'
                 << ecf.p_ << '\t'
                 << ecf.n_ << '\t'
                 << ecf.r_ << '\t'
                 << ecf.eta_;
}

inline
int etafactorvector_toh5file(const std::vector<EtaCreationFactor> efvec,
			     H5::H5File file,
			     std::string sgroup = "Electrostatic_interaction",
			     std::string sdataset = "eta_factor_vector") {

  try
    {    
      // Turn off the auto-printing
      H5::Exception::dontPrint();      
      
      hsize_t dim[] = {efvec.size()};
      H5::DataSpace dspace(1, dim);
      
      H5::CompType etatype(sizeof(EtaCreationFactor));
      
      etatype.insertMember("Id", HOFFSET(EtaCreationFactor, id_), H5::PredType::NATIVE_UINT);
      etatype.insertMember("l", HOFFSET(EtaCreationFactor, l_), H5::PredType::NATIVE_SHORT);
      etatype.insertMember("q", HOFFSET(EtaCreationFactor, q_), H5::PredType::NATIVE_SHORT);
      etatype.insertMember("m", HOFFSET(EtaCreationFactor, m_), H5::PredType::NATIVE_SHORT);
      etatype.insertMember("p", HOFFSET(EtaCreationFactor, p_), H5::PredType::NATIVE_SHORT);
      etatype.insertMember("n", HOFFSET(EtaCreationFactor, n_), H5::PredType::NATIVE_SHORT);
      etatype.insertMember("r", HOFFSET(EtaCreationFactor, r_), H5::PredType::NATIVE_SHORT);
      etatype.insertMember("eta", HOFFSET(EtaCreationFactor, eta_), H5::PredType::NATIVE_DOUBLE);


      H5::Group group;
      // attempt to open group
      try
	{
	  H5::Exception::dontPrint();
	  group = file.openGroup(sgroup);
	}
      catch( H5::Exception gerror ) {
	// if group does not exist, create it
	try
	  {
	    H5::Exception::dontPrint();
	    group = file.createGroup(sgroup);
	  }
	catch( H5::Exception gerror2 ) {
#ifdef H5API110
	  gerror2.printErrorStack();
#else
	  gerror2.printError();
#endif
	  std::cout << std::endl << "[ee] h5 I/O error in group "
		    << sgroup << ". Terminate.\n";
	  std::terminate();
	}
      }
      
      H5::DataSet dataset;
      // attempt to open dataset
      try
	{
	  H5::Exception::dontPrint();
	  dataset = group.openDataSet(sdataset);
	}
      catch( H5::Exception aerror )
	{
	  // create dataset if not exists
	  try
	    {
	      dataset = group.createDataSet(sdataset, etatype, dspace);
	      dspace.close();
	    }
	  catch( H5::Exception aerror2 )
	    {	      
#ifdef H5API110
	      aerror2.printErrorStack();
#else
	      aerror2.printError();
#endif	
	      std::cout << std::endl << "[ee] h5 I/O error in dataset "
			<< sdataset << ". Terminate.\n";
	      std::terminate();
	    }
	}      
      
      // Write data to the dataset, pass plain array
      dataset.write(&efvec[0], etatype);
      
      dataset.close();      
    }// end try block
  // catch failure caused by the H5File operations
  catch( H5::FileIException error )
    {
#ifdef H5API110
      error.printErrorStack();
#else
      error.printError();
#endif
      return -1;
    }
   // catch failure caused by the DataSet operations
  catch( H5::DataSetIException error )
   {
#ifdef H5API110
     error.printErrorStack();
#else
     error.printError();
#endif
     return -1;
   }
  // catch failure caused by the DataSpace operations
  catch( H5::DataSpaceIException error )
    {
#ifdef H5API110
      error.printErrorStack();
#else
      error.printError();
#endif
      return -1;
   }
  // catch failure caused by the DataSpace operations
  catch( H5::DataTypeIException error )
    {
#ifdef H5API110
      error.printErrorStack();
#else
      error.printError();
#endif
      return -1;
    }
  return 0;  
}

inline
void etafactor_fromstream(EtaCreationFactor &ecf, std::istringstream &ss) {
  std::string sid_;           //!< Identifier
  std::string sl_;            //!< Index l of current volume section
  std::string sq_;            //!< Index q of current charge section
  std::string sm_;            //!< Index m of coalescing volume
  std::string sp_;            //!< Index p of coalescing charge
  std::string sn_;            //!< Index n of coalescing volume
  std::string sr_;            //!< Index r of coalescing volume
  std::string seta_;                //!< Value for \f$\eta_{mp,nr}(v_l, Q_q)\f$
  
  unsigned int id_;           //!< Identifier
  short l_;            //!< Index l of current volume section
  short q_;            //!< Index q of current charge section
  short m_;            //!< Index m of coalescing volume
  short p_;            //!< Index p of coalescing charge
  short n_;            //!< Index n of coalescing volume
  short r_;            //!< Index r of coalescing volume
  double eta_;                //!< Value for \f$\eta_{mp,nr}(v_l, Q_q)\f$

  ss >> sid_
     >> sl_
     >> sq_
     >> sm_
     >> sp_
     >> sn_
     >> sr_
     >> seta_;

  id_  = boost::lexical_cast<unsigned int>(sid_ );
  l_   = boost::lexical_cast<short>(sl_  );
  q_   = boost::lexical_cast<short>(sq_  );
  m_   = boost::lexical_cast<short>(sm_  );
  p_   = boost::lexical_cast<short>(sp_  );
  n_   = boost::lexical_cast<short>(sn_  );
  r_   = boost::lexical_cast<short>(sr_  );
  eta_ = std::stod(seta_);

  fill_eta(ecf, id_, l_, q_, m_, p_, n_, r_, eta_);
}


class DeathFactor {
public:

  DeathFactor() {
    // call fill to set 0 all the fields
    fill_death(*this, 0, 0, 0, 0, 0, 0.0);
  }

  unsigned int id_;           //!< Identifier
  short l_;            //!< Index l of current volume section
  short q_;            //!< Index q of current charge section
  short m_;            //!< Index m of coalescing volume
  short p_;            //!< Index p of coalescing charge
  double death_;              //!< Value for death term

  friend void fill_death(DeathFactor &dth,
                       unsigned int id,
                       short l,
                       short q,
                       short m,
                       short p,
                       double death);

  friend void deathfactor_tostream(const DeathFactor &dth,
                                   std::fstream& deathfactor_file);

  friend void deathfactor_fromstream(DeathFactor &df,
                                     std::istringstream &ss);
};

inline
void fill_death(DeathFactor &dth,
                unsigned int id,
                short l,
                short q,
                short m,
                short p,
                double death) {
  dth.id_=id;
  dth.l_=l;
  dth.q_=q;
  dth.m_=m;
  dth.p_=p;
  dth.death_=death;
}

inline
void deathfactor_tostream(const DeathFactor &dth, std::fstream& deathfactor_file) {
  deathfactor_file << dth.id_ << '\t'
                   << dth.l_ << '\t'
                   << dth.q_ << '\t'
                   << dth.m_ << '\t'
                   << dth.p_ << '\t'
                   << dth.death_;
}

inline
void deathfactor_fromstream(DeathFactor &df, std::istringstream &ss) {
  std::string sid_;           //!< Identifier
  std::string sl_;            //!< Index l of current volume section
  std::string sq_;            //!< Index q of current charge section
  std::string sm_;            //!< Index m of coalescing volume
  std::string sp_;            //!< Index p of coalescing charge
  std::string sdeath_;        //!< Value for death term
  
  unsigned int id_;           //!< Identifier
  short l_;            //!< Index l of current volume section
  short q_;            //!< Index q of current charge section
  short m_;            //!< Index m of coalescing volume
  short p_;            //!< Index p of coalescing charge
  double death_;              //!< Value for \f$\eta_{mp,nr}(v_l, Q_q)\f$

  ss >> sid_
     >> sl_
     >> sq_
     >> sm_
     >> sp_
     >> sdeath_;

  id_  = boost::lexical_cast<unsigned int>(sid_ );
  l_   = boost::lexical_cast<short>(sl_  );
  q_   = boost::lexical_cast<short>(sq_  );
  m_   = boost::lexical_cast<short>(sm_  );
  p_   = boost::lexical_cast<short>(sp_  );
  death_ = std::stod(sdeath_);

  fill_death(df, id_, l_, q_, m_, p_, death_);
}


inline
int deathfactorvector_toh5file(const std::vector<DeathFactor> dfvec,
			     H5::H5File file,
			     std::string sgroup = "Electrostatic_interaction",
			     std::string sdataset = "death_factor_vector") {

  try
    {    
      // Turn off the auto-printing
      H5::Exception::dontPrint();      
      
      hsize_t dim[] = {dfvec.size()};
      H5::DataSpace dspace(1, dim);
      
      H5::CompType etatype(sizeof(DeathFactor));
      
      etatype.insertMember("Id", HOFFSET(DeathFactor, id_), H5::PredType::NATIVE_UINT);
      etatype.insertMember("l", HOFFSET(DeathFactor, l_), H5::PredType::NATIVE_SHORT);
      etatype.insertMember("q", HOFFSET(DeathFactor, q_), H5::PredType::NATIVE_SHORT);
      etatype.insertMember("m", HOFFSET(DeathFactor, m_), H5::PredType::NATIVE_SHORT);
      etatype.insertMember("p", HOFFSET(DeathFactor, p_), H5::PredType::NATIVE_SHORT);
      etatype.insertMember("death", HOFFSET(DeathFactor, death_), H5::PredType::NATIVE_DOUBLE);

      H5::Group group;
      // attempt to open group
      try
	{
	  H5::Exception::dontPrint();
	  group = file.openGroup(sgroup);
	}
      catch( H5::Exception gerror ) {
	// if group does not exist, create it
	try
	  {
	    H5::Exception::dontPrint();
	    group = file.createGroup(sgroup);
	  }
	catch( H5::Exception gerror2 ) {
#ifdef H5API110
	  gerror2.printErrorStack();
#else
	  gerror2.printError();
#endif
	  std::cout << std::endl << "[ee] h5 I/O error in group "
		    << sgroup << ". Terminate.\n";
	  std::terminate();
	}
      }
      
      H5::DataSet dataset;
      // attempt to open dataset
      try
	{
	  H5::Exception::dontPrint();
	  dataset = group.openDataSet(sdataset);
	}
      catch( H5::Exception aerror )
	{
	  // create dataset if not exists
	  try
	    {
	      dataset = group.createDataSet(sdataset, etatype, dspace);
	      dspace.close();
	    }
	  catch( H5::Exception aerror2 )
	    {	      
#ifdef H5API110
	      aerror2.printErrorStack();
#else
	      aerror2.printError();
#endif	
	      std::cout << std::endl << "[ee] h5 I/O error in dataset "
			<< sdataset << ". Terminate.\n";
	      std::terminate();
	    }
	}      
      
      // Write data to the dataset, pass plain array
      dataset.write(&dfvec[0], etatype);
      
      dataset.close();      
    }// end try block
  // catch failure caused by the H5File operations
  catch( H5::FileIException error )
    {
#ifdef H5API110
      error.printErrorStack();
#else
      error.printError();
#endif
      return -1;
    }
   // catch failure caused by the DataSet operations
  catch( H5::DataSetIException error )
   {
#ifdef H5API110
     error.printErrorStack();
#else
     error.printError();
#endif
     return -1;
   }
  // catch failure caused by the DataSpace operations
  catch( H5::DataSpaceIException error )
    {
#ifdef H5API110
      error.printErrorStack();
#else
      error.printError();
#endif
      return -1;
   }
  // catch failure caused by the DataSpace operations
  catch( H5::DataTypeIException error )
    {
#ifdef H5API110
      error.printErrorStack();
#else
      error.printError();
#endif
      return -1;
    }
  return 0;  
}
#endif
