/*
 * Copyright 2018 Benjamin Santos <caos21@gmail.com>
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

#ifndef ENHANCEMENT_H
#define ENHANCEMENT_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <utility>
#include <functional>
#include <algorithm>
#include <parallel/algorithm>
#include <cstdlib>
#include <memory>
#include <chrono>
#include <ctime>

#include <omp.h>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include "log.h"
#include "eint.h"
#include "array.h"
#include "constants.h"
#include "utils.h"

#define MPC_ERROR 1.0//0.4

using namespace utils;

namespace enhancement {

  const boost::uintmax_t ROOT_MAXITER = 1000;
  const tools::eps_tolerance<double> ROOT_TOL(30);

  class ParticlePair {
  public:
    //
    ParticlePair() {
      // call fill to set 0 all the fields
      //fill_ppair(*this, 0, 0, 0, 0, 0.0, 0.0, true);
      id_ = 0;
      l_ = 0;
      q_ = 0;
      m_ = 0;
      p_ = 0;
      r21_ = 0.0;
      q21_ = 0.0;
      notswapd_ = true;
    }

    ParticlePair(const unsigned long &id,
		 const short &l,
		 const short &q,
		 const short &m,
		 const short &p,
		 const double &r21,
		 const double &q21,
		 const bool& notswapd) :
      id_(id),l_(l), q_(q), m_(m), p_(p), r21_(r21), q21_(q21), notswapd_(notswapd){
      // call fill to set 0 all the fields
      //fill_ppair(*this, 0, 0, 0, 0, 0.0, 0.0, true);
    }
    
    ParticlePair(const ParticlePair &pp) {
      assign(pp);
    }

    unsigned long id_;   //!< Index
    short l_;            //!< Index l of coalescing particle 1
    short q_;            //!< Index q of coalescing particle 1
    short m_;            //!< Index m of coalescing particle 2
    short p_;            //!< Index p of coalescing particle 2
    double r21_;                //!< Value for radius fraction r2/r1
    double q21_;                //!< Value for radius fraction r2/r1
    bool notswapd_;            //!< false if indices were swapped

    void assign(const ParticlePair &pp) {
      id_=pp.id_;
      l_=pp.l_;
      q_=pp.q_;
      m_=pp.m_;
      p_=pp.p_;
      r21_=pp.r21_;
      q21_=pp.q21_;
      notswapd_ = pp.notswapd_;
    }
  
    friend void fill_ppair(ParticlePair &pp,
			   const unsigned long &id,
			   const short &l,
			   const short &q,
			   const short &m,
			   const short &p,
			   const double &r21,
			   const double &q21,
			   const bool& notswapd);

    friend void particlepair_tostream(const ParticlePair &pp,
				      std::ostream& stream);

    friend void particlepair_fromstream(ParticlePair &pp,
					std::istringstream &stream);
    
    void print(std::ostream& stream = std::cout) {
      particlepair_tostream(*this, stream);
    }
  
  };

  inline
  void fill_ppair(ParticlePair &pp,
		  const unsigned long &id,
		  const short &l,
		  const short &q,
		  const short &m,
		  const short &p,
		  const double &r21,
		  const double &q21,
		  const bool& notswapd) {
    pp.id_=id;
    pp.l_=l;
    pp.q_=q;
    pp.m_=m;
    pp.p_=p;
    pp.r21_=r21;
    pp.q21_=q21;
    pp.notswapd_ = notswapd;
  }

  inline
  void particlepair_tostream(const ParticlePair &pp, std::ostream& stream) {
    stream << pp.id_ << '\t'
	   << pp.l_ << '\t'
	   << pp.q_ << '\t'
	   << pp.m_ << '\t'
	   << pp.p_ << '\t'
	   << pp.r21_ << '\t'
	   << pp.q21_ << '\t'
	   << pp.notswapd_ ;
  }

  inline
  void particlepair_fromstream(ParticlePair &pp,
			       std::istringstream &stream) {
    std::string sid_;            //!< Index l of current volume section
    std::string sl_;            //!< Index l of current volume section
    std::string sq_;            //!< Index q of current charge section
    std::string sm_;            //!< Index m of coalescing volume
    std::string sp_;            //!< Index p of coalescing charge
    std::string sr21_;
    std::string sq21_;
    std::string snotswapd_;

    unsigned long id_;
    short l_;            //!< Index l of coalescing particle 1
    short q_;            //!< Index q of coalescing particle 1
    short m_;            //!< Index m of coalescing particle 2
    short p_;            //!< Index p of coalescing particle 2
    double r21_;                //!< Value for radius fraction r2/r1
    double q21_;                //!< Value for radius fraction r2/r1
    bool notswapd_;            //!< false if indices were swapped

    stream >> sid_
	   >> sl_
	   >> sq_
	   >> sm_
	   >> sp_
	   >> sr21_
	   >> sq21_
	   >> snotswapd_;

    id_  = boost::lexical_cast<unsigned long>(sid_);
    l_   = boost::lexical_cast<short>(sl_);
    q_   = boost::lexical_cast<short>(sq_);
    m_   = boost::lexical_cast<short>(sm_);
    p_   = boost::lexical_cast<short>(sp_);
    r21_ = std::stod(sr21_);
    q21_ = std::stod(sq21_);
    notswapd_  = boost::lexical_cast<bool>(snotswapd_);

    fill_ppair(pp, id_, l_, q_, m_, p_, r21_, q21_, notswapd_);
  }

  class ReducedParticlePair {
  public:
    //
    ReducedParticlePair() {
      // call fill to set 0 all the fields
      fill_rppair(*this, 0, 0.0, 0.0);
    }

    bool operator<(const ReducedParticlePair &rpp) const { 
      if (r21_ != rpp.r21_) {
	return (r21_ < rpp.r21_);
      }
      else {
	return (q21_ < rpp.q21_);
      }
    }    

    inline
    bool operator==(const ReducedParticlePair &rpp) const { 
      if ((r21_==rpp.r21_) && (q21_==rpp.q21_)) {
	return true;
      }
      else {
	return false;
      }
    }

    // friend
    // bool operator == (const ReducedParticlePair& lx,
    // 		    const ReducedParticlePair& rx);

    unsigned long id_;
    double r21_;                //!< Value for radius fraction r2/r1
    double q21_;                //!< Value for radius fraction r2/r1

    std::vector<long> repetitions_;

    ReducedParticlePair(const ReducedParticlePair &rpp) {
      assign(rpp);
    }

    void assign(const ReducedParticlePair &rpp) {
      id_ =rpp.id_;
      r21_=rpp.r21_;
      q21_=rpp.q21_;
      repetitions_  = rpp.repetitions_;
    }
  
    friend void reducedparticlepair_tostream(const ReducedParticlePair &rpp,
					     std::ostream& stream);

    friend void reducedparticlepair_fromstream(ReducedParticlePair &rpp,
					       std::istream& stream);

    void print(std::ostream& stream = std::cout) {
      reducedparticlepair_tostream(*this, stream);
    }  

    friend void fill_rppair(ReducedParticlePair &rpp,
			    unsigned long id_,
			    double r21,
			    double q21);

    friend void fill_rppair(ReducedParticlePair &rpp,
			    unsigned long id_,
			    double r21,
			    double q21,
			    std::vector<long> repetitions);
    
    friend void fill_rppair(ReducedParticlePair &rpp,
			    unsigned long id_,
			    double r21,
			    double q21,
			    long irepeated);

    void fill_repeated(long const &irepeated) {
      repetitions_.push_back(irepeated);
    }
  };

  // bool operator == (const ReducedParticlePair& lx,
  // 		  const ReducedParticlePair& rx) {
  //   // return ((lx.r21_==rx.r21_) && (lx.q21_==rx.q21_));
  //   if ((lx.r21_==rx.r21_) && (lx.q21_==rx.q21_)) {
  //     return true;
  //   }
  //   else {
  //     return false;
  //   }
  // }

  inline
  bool operator == (const ReducedParticlePair& lr,
		    const ParticlePair& rp) {
    // return ((lx.r21_==rx.r21_) && (lx.q21_==rx.q21_));
    if ((lr.r21_==rp.r21_) && (lr.q21_==rp.q21_)) {
      return true;
    }
    else {
      return false;
    }
  }

  inline
  void fill_rppair(ReducedParticlePair &rpp,
		   unsigned long id,
		   double r21,
		   double q21) {
    rpp.id_ =id;
    rpp.r21_=r21;
    rpp.q21_=q21;
  }

  inline
  void fill_rppair(ReducedParticlePair &rpp,
		   unsigned long id,
		   double r21,
		   double q21,
		   std::vector<long> repetitions) {
    rpp.id_ =id;
    rpp.r21_=r21;
    rpp.q21_=q21;
    rpp.repetitions_ = repetitions;    
  }
  
  inline
  void fill_rppair(ReducedParticlePair &rpp,
		   unsigned long id,
		   double r21,
		   double q21,
		   long irepeated) {
    rpp.id_=id;
    rpp.r21_=r21;
    rpp.q21_=q21;
    rpp.repetitions_.push_back(irepeated);
  }

  struct reducedPairsComparison {
    bool operator()(const ReducedParticlePair& lx,
		    const ReducedParticlePair& rx) const {
      if (lx.r21_ != rx.r21_) {
	return (lx.r21_ < rx.r21_);
      }
      return (lx.q21_ < rx.q21_);
    }
  };

  inline
  void reducedparticlepair_tostream(const ReducedParticlePair &rpp,
				    std::ostream& stream) {
    stream << rpp.id_ << '\t' 
	   << rpp.r21_ << '\t'
	   << rpp.q21_ << '\t';
	   // << "repeated = \t";
    for (auto ir: rpp.repetitions_) {
      stream << ir << ",";
    }
  }


  inline
  void reducedparticlepair_fromstream(ReducedParticlePair &rpp,
				      std::istream& stream) {
    std::string sid_;
    std::string sr21_;
    std::string sq21_;

    std::string srep_;

    unsigned long id_;
    double r21_;                //!< Value for radius fraction r2/r1
    double q21_;                //!< Value for radius fraction r2/r1

    std::vector<long> repetitions_;
    
    stream >> sid_
           >> sr21_
	   >> sq21_
	   >> srep_;

    id_  = boost::lexical_cast<unsigned long>(sid_);
    r21_ = std::stod(sr21_);
    q21_ = std::stod(sq21_);

    //split line
    std::vector<std::string> vssrep_;
    boost::split(vssrep_, srep_, boost::is_any_of(","));

    for (std::vector<std::string>::iterator it = vssrep_.begin(); it < vssrep_.end(); ++it) {
      try {
	long trep = (boost::lexical_cast<long>(*it));
	repetitions_.push_back(trep);
      }
      catch(const boost::bad_lexical_cast &) {
	//std::cerr << "\n[ee] Error from lexical cast : " << *it;
	break;
      }
      //std::cout << "\t" << *it;
    }

    fill_rppair(rpp, id_, r21_, q21_, repetitions_);
  }
  
  class ContactPotential {
  public:
    ContactPotential(const unsigned long &id = 0,
		     const double &potential = 0,
		     const short &n = 0) : id_(id),
					   potential_(potential),
					   n_(n) {

    }
  

    friend
    void fill_contact_potential(ContactPotential &cpot,
				const unsigned long &id,
				const double &potential,
				const short &n);

    void print(std::ostream& stream = std::cout) {
      contactpotential_tostream(*this, stream);
    }
  
    friend
    void contactpotential_tostream(const ContactPotential &cpot,
				   std::ostream& stream);
  
    unsigned long id_;
    double potential_;
    short n_;
  };

  inline
  void fill_contact_potential(ContactPotential &cpot,
			      const unsigned long &id,
			      const double &potential,
			      const short &n) {
    cpot.id_ = id;
    cpot.potential_ = potential;
    cpot.n_ = n;
  }

  inline
  void contactpotential_tostream(const ContactPotential &cpot,
				 std::ostream& stream) {
    stream << '\n'
	   << cpot.id_ << '\t'
	   << cpot.potential_ << '\t'
	   << cpot.n_ ;
  }


  struct find_id {
    unsigned long id;
    find_id(const unsigned long & id_) : id(id_) {
    }

    bool operator()(const ContactPotential& cpot) const {
      //return rpp.id_ = id
      //std::vector<long> rptd = rpp.repetitions_;
      //std::vector<long>::iterator it = std::find(rptd.begin(), rptd.end(), id);
      //return std::find(rptd.begin(), rptd.end(), id) != rptd.end();
      return (cpot.id_ == id);
    }
  };

  typedef std::vector<ContactPotential>::iterator ContactPotentialIterator;
  
  class Enhancement {
  public:

    Enhancement(const darray& rarray_, const darray& qarray_,
		boost_array4d_ref efactor_,
		boost_array4d_ref cpotentials_,
		boost_array4d_ref bpotentials_,
		boost_array4d_ref rbarriers_,
		double eps_,
		src::severity_logger< severity_level > lg_) :
      rarray(rarray_), qarray(qarray_), efactor(efactor_),
      cpotentials(cpotentials_), bpotentials(bpotentials_),
      rbarriers(rbarriers_),
      eps(eps_), lg(lg_) {

      rarray_size = static_cast<unsigned short>(rarray.size());
      qarray_size = static_cast<unsigned short>(qarray.size());

      ncombs = ((rarray_size*qarray_size)*(rarray_size*qarray_size+1))/2;
      
      particle_pairs.reserve(ncombs);

      potential_threshold = 0.0;

      nmin = 25;
      nmax = 2000;
      nstep = 5;

      // double kt = Kboltz*temperature;
      // std::cerr << "\n KT "<< kt;
      // std::cerr << '\n';
      // std::cerr << "\n pot barrier "<< exp(-0.0 / kt);
      // double pp = 0.0;
      // std::cerr << "\n pot barrier "<< eta_barrier(-1e-19, pp);
      // std::cerr << '\n';
      // std::cerr << "\n q[1] "<< qarray[1];
      // std::cerr << '\n';
      // std::terminate();
      
    }

    inline
    int write_particlepairs(std::string ppfilename) {
      std::ofstream ppstream(std::string(ppfilename+"_pp.dat"));

      ppstream << "#id\tl\tq\tm\tp\tr21\tq21\tnotswapd";
// #pragma omp parallel for ordered//schedule(nonmonotonic:dynamic)
      for (unsigned int i=0; i<particle_pairs.size(); ++i) {
// #pragma omp ordered
	// {
	  ppstream << '\n';
	  particlepair_tostream(particle_pairs[i], ppstream);
	// }
      }
      ppstream.close();

      // Neutral particles
      std::ofstream npstream(std::string(ppfilename+"_np.dat"));

      npstream << "#id\tl\tq\tm\tp\tr21\tq21\tnotswapd";
// #pragma omp parallel for ordered//schedule(nonmonotonic:dynamic)
      for (unsigned int i=0; i<neutral_pairs.size(); ++i) {
// #pragma omp ordered
	// {
	  npstream << '\n';
	  particlepair_tostream(neutral_pairs[i], npstream);
	// }
      }
      npstream.close();

      // Reduced pairs
      std::ofstream rpstream(std::string(ppfilename+"_rp.dat"));

      rpstream << "#id\tr21\tq21\trepetitions";
// #pragma omp parallel for ordered//schedule(nonmonotonic:dynamic)
      for (unsigned int i=0; i<reduced_pairs.size(); ++i) {
// #pragma omp ordered
	// {
	  rpstream << '\n';
	  reducedparticlepair_tostream(reduced_pairs[i], rpstream);
	// }
      }
      rpstream.close();

      return 0;
    }
    
    inline
    int read_particlepairs(std::string ppfilename) {

      // Read particle pairs
      {
	// Clear vector
	particle_pairs.clear();
	  
	std::ifstream ppstream(std::string(ppfilename+"_pp.dat"));
	
	std::string line;
	unsigned int nlines = 0;
	std::getline(ppstream, line);
	while (std::getline(ppstream, line)) {
	  ParticlePair pp;
	  std::istringstream iss(line);	
	  particlepair_fromstream(pp, iss);
	  particle_pairs.push_back(pp);
	  ++nlines;
	}
	ppstream.close();
	BOOST_LOG_SEV(lg, info) <<  "Particle pair size = " << particle_pairs.size();
      }

      // Read particle pairs
      {
	// Clear vector
	neutral_pairs.clear();
	  
	std::ifstream npstream(std::string(ppfilename+"_np.dat"));

	std::string line;
	unsigned int nlines = 0;
	std::getline(npstream, line);
	while (std::getline(npstream, line)) {
	  ParticlePair np;
	  std::istringstream iss(line);	
	  particlepair_fromstream(np, iss);
	  neutral_pairs.push_back(np);
	  ++nlines;
	}
	npstream.close();
	BOOST_LOG_SEV(lg, info) <<  "Neutral pair size = " << neutral_pairs.size();
      }

      // Read reduced pairs
      {
	// Clear vector
	reduced_pairs.clear();
	  
      	std::ifstream rpstream(std::string(ppfilename+"_rp.dat"));

      	std::string line;
      	unsigned int nlines = 0;
      	std::getline(rpstream, line);
      	while (std::getline(rpstream, line)) {
      	  ReducedParticlePair rp;
      	  std::istringstream iss(line);	
      	  reducedparticlepair_fromstream(rp, iss);
      	  reduced_pairs.push_back(rp);
      	  ++nlines;
      	}
      	rpstream.close();
	BOOST_LOG_SEV(lg, info) <<  "Reduced pair size = " << reduced_pairs.size();
      }

	// for (auto rp: reduced_pairs){
	//   std::cout << std::endl;
	//   rp.print();	  
	// }      
      if (reduced_pairs.size()>0) {
	contact_potentials.resize(reduced_pairs.size());
	attractive_potentials.reserve(reduced_pairs.size());
      }
      
      return 0;
    }
    
    inline
    void compute_reducedpairs() {

      BOOST_LOG_SEV(lg, info) << "Computing particle pairs...";
      auto start = std::chrono::system_clock::now();
      //
      unsigned long idpp = 0;// index for particle pairs
      unsigned long idnp = 0;// index for neutral pairs
      // #pragma omp parallel for collapse(4) schedule(auto)
      for (unsigned int l=0; l<rarray_size; ++l) {
	// iterate in charges particle 1
	for (unsigned int q=0; q<qarray_size; ++q) {
	  // iterate in radii particle 2
	  for (unsigned int m=0; m<rarray_size; ++m) {
	    // iterate in charges particle 2	
	    for (unsigned int p=0; p<qarray_size; ++p) {
	      unsigned int ip1 = q * rarray_size + l;
	      unsigned int ip2 = p * rarray_size + m;
	      // avoid repetitions
	      if((p>=q) && (ip2 >= ip1)) {
		double r1 = rarray[l];  double q1 = qarray[q];
		double r2 = rarray[m];  double q2 = qarray[p];

		unsigned int lp1 = l;  unsigned int qp1 = q;
		unsigned int mp2 = m;  unsigned int pp2 = p;

		bool notswapd = true;

		// swap particles. We want to keep r21 < 1
		if(r2>r1){
		  std::swap(r1, r2);   std::swap(q1, q2);
		  std::swap(lp1, mp2); std::swap(qp1, pp2);
		  notswapd = false;
		}
    
		double r21 = r2/r1;  double q21 = q2/q1;

		// WARNING float comparison
		if(q1==0.0){
		  // q1=0, permute particle 1 with 2
		  std::swap(r1, r2);   std::swap(q1, q2);
		  std::swap(lp1, mp2); std::swap(qp1, pp2);	      
		  q21 = 0.0;
		  notswapd = notswapd | false;
		}

		// if(q2==0.0){
		//   q21 = 0.0;
		// }
	    
		// all neutrals case
		if((q1==0.0)&&(q2==0.0)) {
		  // lp1 = l;  qp1 = q;
		  // mp2 = m;  pp2 = p;
		  // q21 = 0.0;
		  // r21 = rarray[m]/rarray[l];
		  ParticlePair neutralpair(idnp, lp1, qp1, mp2, pp2, r21, q21, notswapd);
		  neutral_pairs.push_back(neutralpair);
		  ++idnp;
		  continue;//break;
		}

		ParticlePair ppair(idpp, lp1, qp1, mp2, pp2, r21, q21, notswapd);
		particle_pairs.push_back(ppair);

		++idpp;
		// write symmetric combination
	      }
	    }
	  }
	}
      }

      auto end = std::chrono::system_clock::now();
	
      std::chrono::duration<double> elapsed_seconds = end-start;
      BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();
  
      BOOST_LOG_SEV(lg, info) << "Pairs size : " << particle_pairs.size();
      BOOST_LOG_SEV(lg, info) << "Neutral pairs size : " << neutral_pairs.size();

      BOOST_LOG_SEV(lg, info) << "Total pairs size : " << particle_pairs.size()+neutral_pairs.size();
      BOOST_LOG_SEV(lg, info) << "Total combinations : " << ncombs;

      BOOST_LOG_SEV(lg, info) << "Computing reduced pairs...";
      
      reduced_pairs.resize(particle_pairs.size());

      start = std::chrono::system_clock::now();
      for (unsigned int i = 0; i<particle_pairs.size(); ++i) {
	reduced_pairs[i].id_  = particle_pairs[i].id_;
	reduced_pairs[i].r21_ = particle_pairs[i].r21_;
	reduced_pairs[i].q21_ = particle_pairs[i].q21_;
      }

      end = std::chrono::system_clock::now();
      elapsed_seconds = end-start;
      BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();


      start = std::chrono::system_clock::now();      
      BOOST_LOG_SEV(lg, info) << "Sorting...";
      //__gnu_parallel::sort(reduced_pairs.begin(), reduced_pairs.end(), reducedPairsComparison());
      std::sort(reduced_pairs.begin(), reduced_pairs.end(), reducedPairsComparison());
      end = std::chrono::system_clock::now();
      elapsed_seconds = end-start;
      BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();


      // ******** new stuff
      // vector sorted with repeated 
      std::vector<ReducedParticlePair> sorted_reducedpairs = reduced_pairs;
      // ******** end new stuff
      
      BOOST_LOG_SEV(lg, info) << "Erasing...";
      reduced_pairs.erase(std::unique(reduced_pairs.begin(), reduced_pairs.end() ), reduced_pairs.end());
      //__gnu_parallel::unique_copy(reduced_pairs.begin(), reduced_pairs.end(), reduced_pairs);
  
      // resize in memory
      std::vector<ReducedParticlePair> tmp = reduced_pairs;
      reduced_pairs.swap(tmp);
      // clear memory for tmp
      std::vector<ReducedParticlePair>().swap(tmp);
    
      //****************************************************************************
      // Add indices of repeated to vector repetitions
      //****************************************************************************

      BOOST_LOG_SEV(lg, info) << "Reduced pairs size : " << reduced_pairs.size();

      start = std::chrono::system_clock::now();
      
      /// new stuff

// #pragma omp parallel for schedule(nonmonotonic:dynamic)// collapse(2)
      long jstart = 0;
      for (unsigned int irp=0; irp<reduced_pairs.size(); ++irp) {
	for (long jpp=jstart; jpp<sorted_reducedpairs.size(); ++jpp) {
	  if(reduced_pairs[irp]==sorted_reducedpairs[jpp]) {
// #pragma omp critical
// 	    {
	      // one thread a at time
	    reduced_pairs[irp].fill_repeated(sorted_reducedpairs[jpp].id_);
	    // }
	  }
	  else {
	    jstart = jpp;
	    break;
	  }
	}
      }

      /// end new stuff

      // this stuff works, slowly
// #pragma omp parallel for schedule(nonmonotonic:dynamic)// collapse(2)
//       for (unsigned int irp=0; irp<reduced_pairs.size(); ++irp) {
// 	for (long jpp=0; jpp<particle_pairs.size(); ++jpp) {
// 	  if(reduced_pairs[irp]==particle_pairs[jpp]) {
// #pragma omp critical
// 	    {
// 	      // one thread a at time
// 	      reduced_pairs[irp].fill_repeated(jpp);
// 	    }
// 	  }
// 	}
//       }
      end = std::chrono::system_clock::now();
      elapsed_seconds = end-start;
      BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();      
 
      //delete(ppair);

      start = std::chrono::system_clock::now();
      // resize contact potentials
      if (reduced_pairs.size()>0) {
	contact_potentials.resize(reduced_pairs.size());
	attractive_potentials.reserve(reduced_pairs.size());
      }

      end = std::chrono::system_clock::now();
      elapsed_seconds = end-start;
      BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();      
      
      BOOST_LOG_SEV(lg, info) << "Done computing pairs.";
      // for (auto rp: reduced_pairs){
      //   std::cout << std::endl;
      //   rp.print();
      // }
    }

    inline
    void compute_coulombpotential_contact() {
      BOOST_LOG_SEV(lg, info)
	<< "Computing Coulomb potential at contact for n pairs : "
	<< reduced_pairs.size();
#pragma omp parallel for ordered//schedule(nonmonotonic:dynamic)
      for (unsigned int i=0; i<reduced_pairs.size(); ++i) {

	unsigned long id = reduced_pairs[i].id_;
	double r21 = reduced_pairs[i].r21_;
	double q21 = reduced_pairs[i].q21_;
 
	double rt = 1.0 + r21;
   
	eint::potential_coulomb_funct pcoulfunct(r21, q21);

	// // ipa at contact 
	double pcoul_rt = pcoulfunct(rt);

#pragma omp critical
	{
	  fill_contact_potential(contact_potentials[i], id, pcoul_rt, 0);
	  //ContactPotential cpot(i, pcoul_rt, 0);
	  //contact_potentials.push_back(cpot);
	  //long ti = reduced_pairs[i].repetitions_[0];
	  //fill_contact_potential(contact_potentials[i], ti, pcoul_rt, 0);
	  // //check attractive potential
	  // if (pcoul_rt<potential_threshold){
	  //   ContactPotential acpot(i, pcoul_rt, 0);
	  //   attractive_potentials.push_back(acpot);
	  // }
	  //std::cout << std::endl << pcoul_rt;
	}
      }

      //attractive_potentials.reserve(1);
      attractive_potentials.resize(0);
      std::vector<ContactPotential> tmp = attractive_potentials;
      attractive_potentials.swap(tmp);
      std::vector<ContactPotential>().swap(tmp);

       
      barrier_potentials.resize(0);
      rbarrier_array.resize(0);

      // BOOST_LOG_SEV(lg, info) << "Attractive potentials size : " << attractive_potentials.size();
      BOOST_LOG_SEV(lg, info) << "Done computing Coulomb potential at contact.";
    }
    
    inline
    void compute_ipapotential_contact() {
      BOOST_LOG_SEV(lg, info)
	<< "Computing ipa potential at contact for n pairs : "
	<< reduced_pairs.size();
// #pragma omp parallel for ordered// schedule(nonmonotonic:dynamic)
      for (unsigned int i=0; i<reduced_pairs.size(); ++i) {

	// id of particle pair (not repeated)
	unsigned long id = reduced_pairs[i].id_;
	double r21 = reduced_pairs[i].r21_;
	double q21 = reduced_pairs[i].q21_;
 
	double rt = 1.0 + r21;
   
	eint::potential_ipa_funct pipafunct(r21, q21, eps);

	// // ipa at contact 
	double pipa_rt = pipafunct(rt);

// #pragma omp ordered// critical
// 	{
	  // contact_potentials has the same index of reduced pairs
	  // and same id of particle pairs
	  fill_contact_potential(contact_potentials[i], id, pipa_rt, 0);
	  //check attractive potential
	  if (pipa_rt<potential_threshold) {
	    // the attractive potential has the id of index of reduced pairs (i)
	    ContactPotential acpot(i, pipa_rt, 0);
	    attractive_potentials.push_back(acpot);
	  }
	// }
      }

      std::vector<ContactPotential> tmp = attractive_potentials;
      attractive_potentials.swap(tmp);
      std::vector<ContactPotential>().swap(tmp);

      barrier_potentials.reserve(attractive_potentials.size());
      rbarrier_array.reserve(attractive_potentials.size());

      BOOST_LOG_SEV(lg, info) << "Attractive potentials size : " << attractive_potentials.size();
      BOOST_LOG_SEV(lg, info) << "Done computing ipa potential at contact.";
    }

    inline
    void compute_ipapotential_barrier() {
      BOOST_LOG_SEV(lg, info) << "Computing ipa potential barrier for n pairs : " << reduced_pairs.size();
// #pragma omp parallel for shared(reduced_pairs) schedule(nonmonotonic:dynamic)
      for (unsigned int i=0; i<attractive_potentials.size(); ++i) {
    
	unsigned int index = attractive_potentials[i].id_;
	//unsigned int index = i;//attractive_potentials[i].id_;
	
	double r21 = reduced_pairs[index].r21_;
	double q21 = reduced_pairs[index].q21_;

	double rmin = 1.0 + r21;
	double rmax = 2.5+r21;

	double min = rmin;
	double max = rmax;
     
	eint::force_ipa_funct fipafunct(r21, q21, eps);
	// force ipa at contact 
	double fipa_rmin = fipafunct(rmin);
	// Force at r max
	double fipa_rmax = fipafunct(rmax);

	// checks if minimum exists
	//	#pragma omp ordered
	if (fipa_rmin*fipa_rmax < 0.0) {
	  // std::cerr << "\n[ii] Mixed phi_rt = " << fipa_rmin << '\t' << fipa_rmax;
	  boost::uintmax_t bmax_iter = ROOT_MAXITER;
	  tools::eps_tolerance<double> tol = ROOT_TOL;

	  std::pair<double, double> pair_fipa;
	  try {
	    pair_fipa = tools::toms748_solve(fipafunct, min, max, fipa_rmin, fipa_rmax, tol, bmax_iter);
	  }
	  catch(const std::exception& exc) {
	    std::cerr << '\n' << exc.what() << '\n';
	    //std::terminate();
	    std::cerr << "\n[ee] No barrier\n";
	  }
	  if(bmax_iter > 990){
	    std::cerr << "\n ERROR max iter " << bmax_iter << "\n\n";
	    std::terminate();
	  }
	  double rbarrier = 0.5*(pair_fipa.first+pair_fipa.second);
	  if(rbarrier>=min){
	    //*********************** USE SAME COEFFICIENTS OF FORCE
	    // ipa functor
	    eint::potential_ipa_funct pipafunct(r21, q21, eps);      
	    // ipa at contact
	    double pipa_rbarrier = pipafunct(rbarrier);
	    // the barrier has the index of the id of reduced pairs
// #pragma omp critical
	    // {
	      ContactPotential barrierpot(index, pipa_rbarrier);
	      barrier_potentials.push_back(barrierpot);
	      rbarrier_array.push_back(rbarrier);
	    // }
	  }
	  else{
	    std::cerr << "\n ERROR Negative rbarrier " << rbarrier << '\n';
	    std::terminate();
	  }	
	} 
      }
 
      std::vector<ContactPotential> tmp = barrier_potentials;
      barrier_potentials.swap(tmp);
      std::vector<ContactPotential>().swap(tmp);

      std::vector<double> tmp2 = rbarrier_array;
      rbarrier_array.swap(tmp2);
      std::vector<double>().swap(tmp2);

      BOOST_LOG_SEV(lg, info) << "Barrier potentials size : " <<barrier_potentials.size();
      BOOST_LOG_SEV(lg, info) << "Done computing ipa potential barrier.";
    }
        
    inline
    void compute_mpcpotential_contact() {
      BOOST_LOG_SEV(lg, info)
	<< "Computing mpc potential at contact for n pairs : "
	<< reduced_pairs.size();      
      // WARNING FIXME
      // std::ofstream outfile("pot1.dat");
      // outfile << "#n\tr21\tq21\tmpc\tpipa\terror_pot\terror\tmsg\n";

#pragma omp parallel for schedule(nonmonotonic:dynamic)
      for (unsigned int i=0; i<reduced_pairs.size(); ++i) {

	unsigned long id = reduced_pairs[i].id_;
	double r21 = reduced_pairs[i].r21_;
	double q21 = reduced_pairs[i].q21_;
 
	double rt = 1.0 + r21;
   
	eint::potential_ipa_funct pipafunct(r21, q21, eps);

	// // ipa at contact 
	double pipa_rt = pipafunct(rt);

	double pcomp = pipa_rt;
	unsigned int initer=0;
      
	for(short n=nmin; n <= nmax; n+=nstep, ++initer){
      
	  // mpc functor
	  eint::potential_mpc_funct pmpcfunct(r21, q21, eps, n);
      
	  // mpc at contact
	  double pmpc_rt = pmpcfunct(rt);

	  double error_comp = max_pct_error(pcomp, pmpc_rt);
	  //double error_pipa = max_pct_error(pipa_rt, pmpc_rt);
	
	  //cerr << "\n[ii] n = " << n << "\t err = " << error_comp;
	  if ((error_comp < MPC_ERROR) && (initer>0)){
// #pragma omp critical
// 	    {
	    fill_contact_potential(contact_potentials[i], id, pmpc_rt, n);
	    // outfile << n
	    // 	      << '\t' << r21 << '\t' << q21
	    // 	      << '\t' << pmpc_rt 
	    // 	      << '\t' << pipa_rt
	    // 	      << '\t' << error_pipa
	    // 	      << '\t' << error_comp << '\n';
	      //check attractive potential
#pragma omp critical
	    {	      
		if (pmpc_rt<potential_threshold){
		  ContactPotential acpot(i, pmpc_rt, n);
		  
		  attractive_potentials.push_back(acpot);
		}
	    
		//break;
		n=nmax;
	    }
	  }
	  else {
	    if (n>nmax-nstep){
	      std::cerr << "\n[ww] Max iterations exceeded\n";
	      fill_contact_potential(contact_potentials[i], i, pmpc_rt, n);
	      // outfile << n
	      // 	<< '\t' << r21 << '\t' << q21
	      // 	<< '\t' << pmpc_rt 
	      // 	<< '\t' << pipa_rt
	      // 	<< '\t' << error_pipa
	      // 	<< '\t' << error_comp
	      // 	<< '\t' << "max_iter\n";
	      //check attractive potential
#pragma omp critical
	      {
		ContactPotential epot(i, pmpc_rt, n);
		error_potentials.push_back(epot);
		if (pmpc_rt<potential_threshold){
		    ContactPotential acpot(i, pmpc_rt, n);
		    attractive_potentials.push_back(acpot);
		}
	      
		//break;
	      }
	    }
	  }
	  pcomp = pmpc_rt;
	}
      }

      //      outfile.close();

      std::vector<ContactPotential> tmp = attractive_potentials;
      attractive_potentials.swap(tmp);
      std::vector<ContactPotential>().swap(tmp);

      barrier_potentials.reserve(attractive_potentials.size());
      rbarrier_array.reserve(attractive_potentials.size());

      BOOST_LOG_SEV(lg, info) << "Attractive potentials size : " << attractive_potentials.size();
      BOOST_LOG_SEV(lg, info) << "Error potentials size : " << error_potentials.size();
      BOOST_LOG_SEV(lg, info) << "Done computing MPC potential at contact.";
    }

    inline
    double eta_attractive(double phimin){
      double kt = Kboltz*temperature;
      return 1.0 - phimin / kt;
    }

    inline
    double eta_repulsive(double phimin){
      double kt = Kboltz*temperature;
      return exp(-phimin / kt);
    }

    inline
    double eta_barrier(double phimin, double phimax){
      double kt = Kboltz*temperature;
      return exp(-phimax / kt) * (1.0 + (phimax-phimin) / kt);
    }

    inline
    void compute_mpcpotential_barrier() {
      BOOST_LOG_SEV(lg, info) << "Computing mpc potential barrier for n pairs : " << attractive_potentials.size();
#pragma omp parallel for shared(reduced_pairs) schedule(nonmonotonic:dynamic)
      for (unsigned int i=0; i<attractive_potentials.size(); ++i) {
    
	unsigned int index = attractive_potentials[i].id_;
	unsigned int nterms = attractive_potentials[i].n_;

	double r21 = reduced_pairs[index].r21_;
	double q21 = reduced_pairs[index].q21_;

	double rmin = 1.0 + r21;
	double rmax = 2.5+r21;//100.0 + r21;

	double min = rmin;
	double max = rmax;
     
	eint::force_mpc_funct fmpcfunct(r21, q21, eps, nterms);
	// force mpc at contact 
	double fmpc_rmin = fmpcfunct(rmin);
	// Force at r max
	double fmpc_rmax = fmpcfunct(rmax);

	// checks if minimum exists
	if (fmpc_rmin*fmpc_rmax < 0.0) {
	  // std::cerr << "\n[ii] Mixed phi_rt = " << fmpc_rmin << '\t' << fmpc_rmax;
	  boost::uintmax_t bmax_iter = ROOT_MAXITER;
	  tools::eps_tolerance<double> tol = ROOT_TOL;

	  std::pair<double, double> pair_fmpc;
	  // try {
	  pair_fmpc = tools::toms748_solve(fmpcfunct, min, max, fmpc_rmin, fmpc_rmax, tol, bmax_iter);
	  //pair_fmpc = tools::bisect(fmpcfunct, min, max, tol, bmax_iter);
	  // }
	  // catch(const std::exception& exc) {
	  //   std::cerr << '\n' << exc.what() << '\n';
	  //   std::terminate();
	  // }
	  if(bmax_iter > 990){
	    std::cerr << "\n ERROR max iter " << bmax_iter << "\n\n";
	    std::terminate();
	  }
	  double rbarrier = 0.5*(pair_fmpc.first+pair_fmpc.second);
	  if(rbarrier>=min){
	    //*********************** USE SAME COEFFICIENTS OF FORCE
#pragma omp critical
	    {
	      // mpc functor
	      eint::potential_mpc_funct pmpcfunct(r21, q21, eps, nterms);      
	      // mpc at contact
	      double pmpc_rbarrier = pmpcfunct(rbarrier);	  
	      ContactPotential barrierpot(index, pmpc_rbarrier, nterms);
	      barrier_potentials.push_back(barrierpot);
	      rbarrier_array.push_back(rbarrier);
	    }
	  }
	  else{
	    std::cerr << "\n ERROR Negative rbarrier " << rbarrier << '\n';
	    std::terminate();
	  }	
	} 
      }
 
      std::vector<ContactPotential> tmp = barrier_potentials;
      barrier_potentials.swap(tmp);
      std::vector<ContactPotential>().swap(tmp);

      std::vector<double> tmp2 = rbarrier_array;
      rbarrier_array.swap(tmp2);
      std::vector<double>().swap(tmp2);

      BOOST_LOG_SEV(lg, info) << "Barrier potentials size : " <<barrier_potentials.size();      
      BOOST_LOG_SEV(lg, info) << "Done computing mpc potential barrier.";
    }

    inline
    void compute_enhancement_factor() {

      for (unsigned int l=0; l<rarray_size; ++l) {
	// iterate in charges particle 1
	for (unsigned int q=0; q<qarray_size; ++q) {
	  // iterate in radii particle 2
	  for (unsigned int m=0; m<rarray_size; ++m) {
	    // iterate in charges particle 2	
	    for (unsigned int p=0; p<qarray_size; ++p) {
	      efactor[l][q][m][p] = 0.0;
	      cpotentials[l][q][m][p] = 0.0;
	      bpotentials[l][q][m][p] = 0.0;	      
	    }
	  }
	}
      }
    
      BOOST_LOG_SEV(lg, info) << "Size for array efactor : " << efactor.size();
            
      BOOST_LOG_SEV(lg, info) << "Computing Enhancement Factor";

      // iterate in reduced_pairs
      for (unsigned int irp=0; irp<reduced_pairs.size(); ++irp) {
	
	// this is the index of reduced pairs, i.e., non repeated pair
	unsigned int index = reduced_pairs[irp].id_;

	// contact potential for this reduced pairhas the same ordering
	// of reduced pairs
	double contact_potential = contact_potentials[irp].potential_;
	
	// get repeated combinations (of particle_pairs)
	std::vector<long> reps = reduced_pairs[irp].repetitions_;

	// find if index is in barrier_potentials, i.e., if a potential
	// barrier exists
	ContactPotentialIterator barrier_it = std::find_if(barrier_potentials.begin(),
							   barrier_potentials.end(),
							   find_id(irp));
	
        // iterate in vector of repeated combinations
	for(unsigned int jrep=0; jrep<reps.size(); ++jrep) {
	  long rep_index = reps[jrep];
	  short l = particle_pairs[rep_index].l_;
	  short q = particle_pairs[rep_index].q_;
	  short m = particle_pairs[rep_index].m_;
	  short p = particle_pairs[rep_index].p_;

	  double tr1 = rarray[l];//(!notswapd ? rarray[l] : rarray[m]);
	  //std::cerr << "\nn q is " << q;
	  double tq1 = qarray[q];//(!notswapd ? qarray[q] : qarray[p]);
	  double potprefactor = tq1*tq1/tr1;
	  double phimin = potprefactor*contact_potential;
	  //std::cerr << "\n\n";

	  // we have potential barrier
	  // ** if barrier compute eta       
	  if(barrier_it != barrier_potentials.end()) {
	    // std::cerr << "\n*********************";
	    // std::cerr << "\n* " << (*it).id_;
	    // 			   std::cerr << "\n* " << (*it).potential_;
	    //
	    double phimax = potprefactor*(*barrier_it).potential_;
	    // std::cerr << "\n*********************" << rep_index;
	    // std::cerr << "\n(" << l << ", " << q << ") + (" << m << ", " << p << ")";
	    // std::cerr << "\n*******contact**********\t" << contact_potential;
	    // std::cerr << "\n*******barrier**********\t" << (*barrier_it).potential_;  
	    // std::cerr << "\n*******phimin**********\t" << phimin;
	    // std::cerr << "\n*******phimax**********\t" << phimax;
	    // double eta = eta_barrier(phimin, phimax);
	    // std::cerr << "\n*******eta**********\t" << eta
	    // 	      << '\t' << tr1
	    // 	      << '\t' << tq1
	    // 	      << '\t' << rarray[m]
	    // 	      << '\t' << qarray[p]
	    // 	      << '\t' << particle_pairs[rep_index].notswapd_;
	    // #pragma omp atomic write
	    // efactor[l][q][m][p] = eta;
	    // #pragma omp atomic write
	    // efactor[m][p][l][q] = eta;
	    // update potentials
	    cpotentials[m][p][l][q] = phimin;
	    cpotentials[l][q][m][p] = phimin;
	    bpotentials[m][p][l][q] = phimax;
	    bpotentials[l][q][m][p] = phimax;
	    unsigned int idr = barrier_it - barrier_potentials.begin();
	    rbarriers[l][q][m][p] = rbarrier_array[idr]*tr1;
	    rbarriers[m][p][l][q] = rbarrier_array[idr]*tr1;
	  }
	  else {
	    // ** if not barrier
	    // ** attraction and repulsion
	    if(contact_potential <= potential_threshold){
	      //std::cerr << "\n--------------------" << rep_index;
	      //std::cerr << "\n(" << l << ", " << q << ") + (" << m << ", " << p << ")";
	      //double phimin = qarray[]
	      //std::cerr << "\n----------contact------\t" << contact_potential;
	      //std::cerr << "\n----------phimin----------\t" << phimin;
	      // double eta = eta_attractive(phimin);
	      //std::cerr << "\n----------eta------\t" << eta;
	      // #pragma omp atomic write
	      // efactor[l][q][m][p] = eta;
	      // #pragma omp atomic write
	      // efactor[m][p][l][q] = eta;
	      // update potentials
	      cpotentials[m][p][l][q] = phimin;
	      cpotentials[l][q][m][p] = phimin;
	      //bpotentials[m][p][l][q] = 0.0;
	      //bpotentials[l][q][m][p] = 0.0;
	    }
	    else {
	      // if(contact_potential > potential_threshold ){
	      	//std::cerr << "\n++++++++++++++++++++" << rep_index;
	      	//std::cerr << "\n(" << l << ", " << q << ") + (" << m << ", " << p << ")";
	      	//std::cerr << "\n++++++++++contact++++++++++\t" << contact_potential;
	      	//std::cerr << "\n++++++++++phimin++++++++++\t" << phimin;
	      	// double eta = eta_repulsive(phimin);
	      	//std::cerr << "\n++++++++++eta++++++++++\t" << eta;
	      	// #pragma omp atomic write
	      	// efactor[l][q][m][p] = eta;
	      	// #pragma omp atomic write
	      	//efactor[m][p][l][q] = eta;
		// update potentials
		cpotentials[m][p][l][q] = phimin;
		cpotentials[l][q][m][p] = phimin;
		//bpotentials[m][p][l][q] = 0.0;
		//bpotentials[l][q][m][p] = 0.0;		
	      // }
	    }
	  }

	  //std::cout << std::endl << l << '\t' << q << '\t' << m << '\t' << p << '\t' << efactor[l][q][m][p];
	 
	}
      }

      //#pragma omp parallel for 
      for (unsigned int l=0; l<rarray_size; ++l) {
	// iterate in charges particle 1
	for (unsigned int q=0; q<qarray_size; ++q) {
	  // iterate in radii particle 2
	  for (unsigned int m=0; m<rarray_size; ++m) {
	    // iterate in charges particle 2	
	    for (unsigned int p=0; p<qarray_size; ++p) {
	      double phimin = cpotentials[l][q][m][p];
	      if (phimin <= potential_threshold) {
		double phimax = bpotentials[l][q][m][p];
		double eta = eta_barrier(phimin, phimax);
                //#pragma omp atomic write
		efactor[l][q][m][p] = eta;		
	      }
	      else {
		double eta = eta_repulsive(phimin);
		//#pragma omp atomic write
		efactor[l][q][m][p] = eta;
	      }
	    }
	  }
	}
      }

      BOOST_LOG_SEV(lg, info) << "Computing neutral pairs : " << neutral_pairs.size();
      
      auto start = std::chrono::system_clock::now();
      for (unsigned int i = 0; i<neutral_pairs.size(); ++i) {
      	short l = neutral_pairs[i].l_;
      	short q = neutral_pairs[i].q_;
      	short m = neutral_pairs[i].m_;
      	short p = neutral_pairs[i].p_;
      	efactor[l][q][m][p] = 1.0;
      	efactor[m][p][l][q] = 1.0;
      }

      auto end = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds = end-start;
      BOOST_LOG_SEV(lg, info) << "Done computing neutral pairs...";
      BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();   
      
    }
    //
    darray rarray;    
    darray qarray;

    boost_array4d_ref efactor;
    boost_array4d_ref cpotentials;
    boost_array4d_ref bpotentials;
    boost_array4d_ref rbarriers;
    
    double eps;

    src::severity_logger< severity_level > lg;
    
    unsigned short rarray_size;    
    unsigned short qarray_size;
    unsigned long ncombs;

    short nmin;
    short nmax;
    short nstep;
    
    std::vector<ParticlePair> particle_pairs;
    std::vector<ParticlePair> neutral_pairs;    
    std::vector<ReducedParticlePair> reduced_pairs;
    std::vector<ContactPotential> contact_potentials;
    std::vector<ContactPotential> attractive_potentials;
    std::vector<ContactPotential> error_potentials;  
    std::vector<ContactPotential> barrier_potentials;
    
    double potential_threshold;
    double temperature = 300.0;
    std::vector<double> rbarrier_array;
  };
}

#endif//ENHANCEMENT_H
