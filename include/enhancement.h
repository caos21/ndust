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
#include <utility>
#include <functional>
#include <algorithm>
#include <parallel/algorithm>
#include <cstdlib>

#include <omp.h>

#include "eint.h"
#include "array.h"
#include "constants.h"
#include "utils.h"

using namespace utils;

namespace enhancement {

  class ParticlePair {
  public:
    //
    ParticlePair() {
      // call fill to set 0 all the fields
      fill_ppair(*this, 0, 0, 0, 0, 0.0, 0.0);
    }

    ParticlePair(const ParticlePair &pp) {
      assign(pp);
    }
  
    short l_;            //!< Index l of coalescing particle 1
    short q_;            //!< Index q of coalescing particle 1
    short m_;            //!< Index m of coalescing particle 2
    short p_;            //!< Index p of coalescing particle 2
    double r21_;                //!< Value for radius fraction r2/r1
    double q21_;                //!< Value for radius fraction r2/r1

    void assign(const ParticlePair &pp) {
      l_=pp.l_;
      q_=pp.q_;
      m_=pp.m_;
      p_=pp.p_;
      r21_=pp.r21_;
      q21_=pp.q21_;
    }
  
    friend void fill_ppair(ParticlePair &pp,
			   const short &l,
			   const short &q,
			   const short &m,
			   const short &p,
			   const double &r21,
			   const double &q21);

    friend void particlepair_tostream(const ParticlePair &pp,
				      std::ostream& stream);

    void print(std::ostream& stream = std::cout) {
      particlepair_tostream(*this, stream);
    }
  
  };

  inline
  void fill_ppair(ParticlePair &pp,
		  const short &l,
		  const short &q,
		  const short &m,
		  const short &p,
		  const double &r21,
		  const double &q21) {
    pp.l_=l;
    pp.q_=q;
    pp.m_=m;
    pp.p_=p;
    pp.r21_=r21;
    pp.q21_=q21;
  }

  inline
  void particlepair_tostream(const ParticlePair &pp, std::ostream& stream) {
    stream << '\n'
	   << pp.l_ << '\t'
	   << pp.q_ << '\t'
	   << pp.m_ << '\t'
	   << pp.p_ << '\t'
	   << pp.r21_ << '\t'
	   << pp.q21_;
  }


  class ReducedParticlePair {
  public:
    //
    ReducedParticlePair() {
      // call fill to set 0 all the fields
      fill_rppair(*this, 0.0, 0.0);
    }

    bool operator<(const ReducedParticlePair &rpp) const { 
      if (r21_ != rpp.r21_) {
	return (r21_ < rpp.r21_);
      }
      return (q21_ < rpp.q21_);
    }    

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
    double r21_;                //!< Value for radius fraction r2/r1
    double q21_;                //!< Value for radius fraction r2/r1

    std::vector<long> repetitions_;

    ReducedParticlePair(const ReducedParticlePair &rpp) {
      assign(rpp);
    }

    void assign(const ReducedParticlePair &rpp) {
      r21_=rpp.r21_;
      q21_=rpp.q21_;
      repetitions_  = rpp.repetitions_;
    }
  
    friend void reducedparticlepair_tostream(const ReducedParticlePair &rpp,
					     std::ostream& stream);

    void print(std::ostream& stream = std::cout) {
      reducedparticlepair_tostream(*this, stream);
    }  

    friend void fill_rppair(ReducedParticlePair &rpp,
			    double r21,
			    double q21);

    friend void fill_rppair(ReducedParticlePair &rpp,
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
		   double r21,
		   double q21) {
    rpp.r21_=r21;
    rpp.q21_=q21;
  }

  inline
  void fill_rppair(ReducedParticlePair &rpp,
		   double r21,
		   double q21,
		   long irepeated) {
    rpp.r21_=r21;
    rpp.q21_=q21;
    rpp.repetitions_.push_back(irepeated);
  }

  struct reducedPairsComparison {
    bool operator()(const ReducedParticlePair& lx, const ReducedParticlePair& rx) const {
      if (lx.r21_ != rx.r21_) {
	return (lx.r21_ < rx.r21_);
      }
      return (lx.q21_ < rx.q21_);
    }
  };

  inline
  void reducedparticlepair_tostream(const ReducedParticlePair &rpp,
				    std::ostream& stream) {
    stream << '\n'
	   << rpp.r21_ << '\t'
	   << rpp.q21_ << '\t'
	   << "repeated = \t";
    for (auto ir: rpp.repetitions_) {
      stream << ir << ", ";
    }
  }

  class ContactPotential {
  public:
    ContactPotential(const long &id = 0,
		     const double &potential = 0,
		     const short &n = 0) : id_(id),
					   potential_(potential),
					   n_(n) {
    }
  

    friend
    void fill_contact_potential(ContactPotential &cpot,
				const long &id,
				const double &potential,
				const short &n);

    void print(std::ostream& stream = std::cout) {
      contactpotential_tostream(*this, stream);
    }
  
    friend
    void contactpotential_tostream(const ContactPotential &cpot,
				   std::ostream& stream);
  
    long id_;
    double potential_;
    long n_;
  };

  inline
  void fill_contact_potential(ContactPotential &cpot,
			      const long &id,
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


  class Enhancement {
  public:

    Enhancement(const darray& rarray_, const darray& qarray_, double eps_) :
      rarray(rarray_), qarray(qarray_), eps(eps_) {

      rarray_size = rarray.size();
      qarray_size = qarray.size();

      ncombs = ((rarray_size*qarray_size)*(rarray_size*qarray_size+1))/2;
      
      particle_pairs.reserve(ncombs);

      potential_threshold = 0.0;
    }
    
    inline
    void compute_reducedpairs() {

      ParticlePair *ppair = new ParticlePair();

      unsigned int i = 0;
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
	    
		if(r2>r1){
		  std::swap(r1, r2);   std::swap(q1, q2);
		  std::swap(lp1, mp2); std::swap(qp1, pp2);
		}
    
		double r21 = r2/r1;  double q21 = q2/q1;

		// WARNING float comparison
		if(q1==0.0){
		  // q1=0, permute particle 1 with 2
		  std::swap(r1, r2);   std::swap(q1, q2);
		  std::swap(lp1, mp2); std::swap(qp1, pp2);	      
		  q21 = 0.0;
		}

		// if(q2==0.0){
		//   q21 = 0.0;
		// }
	    
		// all neutrals case
		if((q1==0.0)&&(q2==0.0)) {
		  lp1 = l;  qp1 = q;
		  mp2 = m;  pp2 = p;
		  q21 = 0.0;
		  r21 = rarray[m]/rarray[l];
		}

		fill_ppair(*ppair, lp1, qp1, mp2, pp2, r21, q21);
		particle_pairs.push_back(*ppair);
	    
		++i;
		// write symmetric combination
	      }
	    }
	  }
	}
      }

      // fill 
      // ReducedParticlePair *rppair = new ReducedParticlePair();
      // fill_rppair(*rppair, particle_pairs[0].r21_, particle_pairs[0].q21_, 0); 
      // reduced_pairs.push_back(*rppair);

      reduced_pairs.resize(particle_pairs.size());

      for (unsigned int i = 0; i<particle_pairs.size(); ++i) {
	reduced_pairs[i].r21_ = particle_pairs[i].r21_;
	reduced_pairs[i].q21_ = particle_pairs[i].q21_;
      }

      //  __gnu_parallel::sort(reduced_pairs.begin(), reduced_pairs.end(), reducedPairsComparison());
      std::sort(reduced_pairs.begin(), reduced_pairs.end(), reducedPairsComparison());

      reduced_pairs.erase(std::unique(reduced_pairs.begin(), reduced_pairs.end() ), reduced_pairs.end());
      //__gnu_parallel::unique_copy(reduced_pairs.begin(), reduced_pairs.end(), reduced_pairs);
  
      // resize in memory
      std::vector<ReducedParticlePair> tmp = reduced_pairs;
      reduced_pairs.swap(tmp);
      // clear memory for tmp
      std::vector<ReducedParticlePair>().swap(tmp);
    
      //****************************************************************************
      // Add indices of repeated to vector ofast repetitions
      //****************************************************************************


#pragma omp parallel for schedule(nonmonotonic:dynamic)
      for (unsigned int irp=0; irp<reduced_pairs.capacity(); ++irp) {
	for (unsigned int jpp=0; jpp<particle_pairs.capacity(); ++jpp) {
	  // #pragma omp critical
	  // {
	  if(reduced_pairs[irp]==particle_pairs[jpp]) {
	    //std::cout << "\n[ii] equals: " << irp << '\t' << jpp;
#pragma omp critical
	    {	    
	      reduced_pairs[irp].fill_repeated(jpp);
	    }
	  }
	  // }
	}
      }
 
      //delete(rppair);
      delete(ppair);

      // resize contact potentials
      contact_potentials.resize(reduced_pairs.size());
      attractive_potentials.reserve(reduced_pairs.size());
      
    }

    inline
    void compute_contact_potentials() {
      std::ofstream outfile("pot1.dat");
      outfile << "#n\tr21\tq21\tmpc\tpipa\terror_pot\terror\tmsg\n";
 
      unsigned int nmin = 25;//25;//25;//10;//80;//20;
      unsigned int nmax = 300;//100;//100;//90;//90;//32;
      unsigned int nstep = 5;//20;

      unsigned int i;
#pragma omp parallel for private(i) schedule(nonmonotonic:dynamic)
      for (i=0; i<reduced_pairs.size(); ++i) {


	double r21 = reduced_pairs[i].r21_;
	double q21 = reduced_pairs[i].q21_;

	// for(auto ups: reduced_pairs){

	//   double r21 = ups.r21_;
	//   double q21 = ups.q21_;
 
	double rt = 1.0 + r21;
   
	eint::potential_ipa_funct pipafunct(r21, q21, eps);

	// // ipa at contact 
	double pipa_rt = pipafunct(rt);

	double pcomp = pipa_rt;
	unsigned int initer=0;
      
	for(unsigned int n=nmin; n <= nmax; n+=nstep, ++initer){
      
	  // mpc functor
	  eint::potential_mpc_funct pmpcfunct(r21, q21, eps, n);
      
	  // mpc at contact
	  double pmpc_rt = pmpcfunct(rt);

	  double error_comp = max_pct_error(pcomp, pmpc_rt);
	  double error_pipa = max_pct_error(pipa_rt, pmpc_rt);
	
	  //cerr << "\n[ii] n = " << n << "\t err = " << error_comp;
	  if ((error_comp < 2.0) && (initer>0)){
#pragma omp critical
	    {
	      fill_contact_potential(contact_potentials[i], i, pmpc_rt, n);
	      outfile << n
		      << '\t' << r21 << '\t' << q21
		      << '\t' << pmpc_rt 
		      << '\t' << pipa_rt
		      << '\t' << error_pipa
		      << '\t' << error_comp << '\n';
	      //check attractive potential
	      if (pmpc_rt<potential_threshold){
		ContactPotential acpot(i, pmpc_rt, n);
		attractive_potentials.emplace_back(acpot);
	      }
	    }
	    break;
	  }
	  else {
	    if (n>nmax-nstep){
#pragma omp critical
	      {
		std::cerr << "\n[ww] Max iterations exceeded\n";
		fill_contact_potential(contact_potentials[i], i, pmpc_rt, n);
		outfile << n
			<< '\t' << r21 << '\t' << q21
			<< '\t' << pmpc_rt 
			<< '\t' << pipa_rt
			<< '\t' << error_pipa
			<< '\t' << error_comp
			<< '\t' << "max_iter\n";
		//check attractive potential
		if (pmpc_rt<potential_threshold){
		  ContactPotential acpot(i, pmpc_rt, n);
		  attractive_potentials.push_back(acpot);
		}
	      }
	      break;
	    }
	  }
	  pcomp = pmpc_rt;
	}
      }

      outfile.close();

  
      std::cout << "\nPotentials attcontact contact : " << attractive_potentials.size() <<  '\n';
      std::vector<ContactPotential> tmp = attractive_potentials;
      std::cout << "\nPotentials attcontact capacity : " << attractive_potentials.capacity() <<  '\n';
      attractive_potentials.swap(tmp);
      std::cout << "\nPotentials attcontact capacity : " << attractive_potentials.capacity() <<  '\n';
      std::vector<ContactPotential>().swap(tmp);
      
    }
    darray rarray;
    
    darray qarray;

    double eps;

    unsigned short rarray_size;
    
    unsigned short qarray_size;

    unsigned long ncombs;
    
    std::vector<ParticlePair> particle_pairs;

    std::vector<ReducedParticlePair> reduced_pairs;

    std::vector<ContactPotential> contact_potentials;

    std::vector<ContactPotential> attractive_potentials;

    double potential_threshold;
  };
}

#endif//ENHANCEMENT_H
