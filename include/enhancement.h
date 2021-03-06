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

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <parallel/algorithm>
#include <vector>
#include <sstream>
#include <string>
#include <utility>

#include <omp.h>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "array.h"
#include "constants.h"
#include "eint.h"
#include "log.h"
#include "utils.h"

#define MPC_ERROR 1.0  // 0.4

#define ETA_TOL 1.0e-16

using namespace utils;


typedef std::pair<short, short> pp;

typedef std::pair<pp, pp> pairs;

typedef std::vector<short> vs;

typedef std::vector<pp> vp;

typedef std::vector<pairs> vpairs;

namespace enhancement {

const boost::uintmax_t ROOT_MAXITER = 1000;
const tools::eps_tolerance<double> ROOT_TOL(30);

/*!
 * class PairElement
 *
 * Represents a value for a particle pair combination (l,q)+(m,p)
 *
 */
class PairElement {
 public:
  //
  PairElement() {
    l_ = 0;
    q_ = 0;
    m_ = 0;
    p_ = 0;
    value_ = 0.0;
  }

  PairElement(const short &l, const short &q, const short &m, const short &p,
              const double &value)
      : l_(l), q_(q), m_(m), p_(p), value_(value) {}

  PairElement(const PairElement &pe) { assign(pe); }

  short l_;       //!< Index l of coalescing particle 1
  short q_;       //!< Index q of coalescing particle 1
  short m_;       //!< Index m of coalescing particle 2
  short p_;       //!< Index p of coalescing particle 2
  double value_;  //!< Value

  void assign(const PairElement &pe) {
    l_ = pe.l_;
    q_ = pe.q_;
    m_ = pe.m_;
    p_ = pe.p_;
    value_ = pe.value_;
  }

  friend void fill_pe(PairElement &pe, const short &l, const short &q,
                      const short &m, const short &p, const double &value);

  friend void pairelement_tostream(const PairElement &pe, std::ostream &stream);

  friend void pairelement_tostream(PairElement &pe, std::istringstream &stream);

  void print(std::ostream &stream = std::cout) {
    pairelement_tostream(*this, stream);
  }
};

inline void get_pe(const PairElement &pe, short &l, short &q, short &m,
                   short &p, double &value) {
  l = pe.l_;
  q = pe.q_;
  m = pe.m_;
  p = pe.p_;
  value = pe.value_;
}

inline void fill_pe(PairElement &pe, const short &l, const short &q,
                    const short &m, const short &p, const double &value) {
  pe.l_ = l;
  pe.q_ = q;
  pe.m_ = m;
  pe.p_ = p;
  pe.value_ = value;
}

inline void pairelement_tostream(const PairElement &pe, std::ostream &stream) {
  stream << pe.l_ << '\t' << pe.q_ << '\t' << pe.m_ << '\t' << pe.p_ << '\t'
         << pe.value_;
}

inline void pairelement_fromstream(PairElement &pe,
                                   std::istringstream &stream) {
  std::string sl_;  //!< Index l of current volume section
  std::string sq_;  //!< Index q of current charge section
  std::string sm_;  //!< Index m of coalescing volume
  std::string sp_;  //!< Index p of coalescing charge
  std::string svalue_;

  short l_;       //!< Index l of coalescing particle 1
  short q_;       //!< Index q of coalescing particle 1
  short m_;       //!< Index m of coalescing particle 2
  short p_;       //!< Index p of coalescing particle 2
  double value_;  //!< Value

  stream >> sl_ >> sq_ >> sm_ >> sp_ >> svalue_;

  l_ = boost::lexical_cast<short>(sl_);
  q_ = boost::lexical_cast<short>(sq_);
  m_ = boost::lexical_cast<short>(sm_);
  p_ = boost::lexical_cast<short>(sp_);
  value_ = std::stod(svalue_);

  fill_pe(pe, l_, q_, m_, p_, value_);
}

//

class ParticlePair {
 public:
  //
  ParticlePair() {
    // call fill to set 0 all the fields
    // fill_ppair(*this, 0, 0, 0, 0, 0.0, 0.0, true);
    id_ = 0;
    l_ = 0;
    q_ = 0;
    m_ = 0;
    p_ = 0;
    r21_ = 0.0;
    q21_ = 0.0;
    notswapd_ = true;
  }

  ParticlePair(const unsigned long &id, const short &l, const short &q,
               const short &m, const short &p, const double &r21,
               const double &q21, const bool &notswapd)
      : id_(id),
        l_(l),
        q_(q),
        m_(m),
        p_(p),
        r21_(r21),
        q21_(q21),
        notswapd_(notswapd) {
    // call fill to set 0 all the fields
    // fill_ppair(*this, 0, 0, 0, 0, 0.0, 0.0, true);
  }

  ParticlePair(const ParticlePair &pp) { assign(pp); }

  unsigned long id_;  //!< Index
  short l_;           //!< Index l of coalescing particle 1
  short q_;           //!< Index q of coalescing particle 1
  short m_;           //!< Index m of coalescing particle 2
  short p_;           //!< Index p of coalescing particle 2
  double r21_;        //!< Value for radius fraction r2/r1
  double q21_;        //!< Value for radius fraction r2/r1
  bool notswapd_;     //!< false if indices were swapped

  void assign(const ParticlePair &pp) {
    id_ = pp.id_;
    l_ = pp.l_;
    q_ = pp.q_;
    m_ = pp.m_;
    p_ = pp.p_;
    r21_ = pp.r21_;
    q21_ = pp.q21_;
    notswapd_ = pp.notswapd_;
  }

  friend void fill_ppair(ParticlePair &pp, const unsigned long &id,
                         const short &l, const short &q, const short &m,
                         const short &p, const double &r21, const double &q21,
                         const bool &notswapd);

  friend void particlepair_tostream(const ParticlePair &pp,
                                    std::ostream &stream);

  friend void particlepair_fromstream(ParticlePair &pp,
                                      std::istringstream &stream);

  void print(std::ostream &stream = std::cout) {
    particlepair_tostream(*this, stream);
  }
};

inline void fill_ppair(ParticlePair &pp, const unsigned long &id,
                       const short &l, const short &q, const short &m,
                       const short &p, const double &r21, const double &q21,
                       const bool &notswapd) {
  pp.id_ = id;
  pp.l_ = l;
  pp.q_ = q;
  pp.m_ = m;
  pp.p_ = p;
  pp.r21_ = r21;
  pp.q21_ = q21;
  pp.notswapd_ = notswapd;
}

inline void particlepair_tostream(const ParticlePair &pp,
                                  std::ostream &stream) {
  stream << pp.id_ << '\t' << pp.l_ << '\t' << pp.q_ << '\t' << pp.m_ << '\t'
         << pp.p_ << '\t' << pp.r21_ << '\t' << pp.q21_ << '\t' << pp.notswapd_;
}

inline void particlepair_fromstream(ParticlePair &pp,
                                    std::istringstream &stream) {
  std::string sid_;  //!< Index l of current volume section
  std::string sl_;   //!< Index l of current volume section
  std::string sq_;   //!< Index q of current charge section
  std::string sm_;   //!< Index m of coalescing volume
  std::string sp_;   //!< Index p of coalescing charge
  std::string sr21_;
  std::string sq21_;
  std::string snotswapd_;

  unsigned long id_;
  short l_;        //!< Index l of coalescing particle 1
  short q_;        //!< Index q of coalescing particle 1
  short m_;        //!< Index m of coalescing particle 2
  short p_;        //!< Index p of coalescing particle 2
  double r21_;     //!< Value for radius fraction r2/r1
  double q21_;     //!< Value for radius fraction r2/r1
  bool notswapd_;  //!< false if indices were swapped

  stream >> sid_ >> sl_ >> sq_ >> sm_ >> sp_ >> sr21_ >> sq21_ >> snotswapd_;

  id_ = boost::lexical_cast<unsigned long>(sid_);
  l_ = boost::lexical_cast<short>(sl_);
  q_ = boost::lexical_cast<short>(sq_);
  m_ = boost::lexical_cast<short>(sm_);
  p_ = boost::lexical_cast<short>(sp_);
  r21_ = std::stod(sr21_);
  q21_ = std::stod(sq21_);
  notswapd_ = boost::lexical_cast<bool>(snotswapd_);

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
    } else {
      return (q21_ < rpp.q21_);
    }
  }

  inline bool operator==(const ReducedParticlePair &rpp) const {
    if ((r21_ == rpp.r21_) && (q21_ == rpp.q21_)) {
      return true;
    } else {
      return false;
    }
  }

  // friend
  // bool operator == (const ReducedParticlePair& lx,
  // 		    const ReducedParticlePair& rx);

  unsigned long id_;
  double r21_;  //!< Value for radius fraction r2/r1
  double q21_;  //!< Value for radius fraction r2/r1

  std::vector<long> repetitions_;

  ReducedParticlePair(const ReducedParticlePair &rpp) { assign(rpp); }

  void assign(const ReducedParticlePair &rpp) {
    id_ = rpp.id_;
    r21_ = rpp.r21_;
    q21_ = rpp.q21_;
    repetitions_ = rpp.repetitions_;
  }

  friend void reducedparticlepair_tostream(const ReducedParticlePair &rpp,
                                           std::ostream &stream);

  friend void reducedparticlepair_fromstream(ReducedParticlePair &rpp,
                                             std::istream &stream);

  void print(std::ostream &stream = std::cout) {
    reducedparticlepair_tostream(*this, stream);
  }

  friend void fill_rppair(ReducedParticlePair &rpp, unsigned long id_,
                          double r21, double q21);

  friend void fill_rppair(ReducedParticlePair &rpp, unsigned long id_,
                          double r21, double q21,
                          std::vector<long> repetitions);

  friend void fill_rppair(ReducedParticlePair &rpp, unsigned long id_,
                          double r21, double q21, long irepeated);

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

inline bool operator==(const ReducedParticlePair &lr, const ParticlePair &rp) {
  // return ((lx.r21_==rx.r21_) && (lx.q21_==rx.q21_));
  if ((lr.r21_ == rp.r21_) && (lr.q21_ == rp.q21_)) {
    return true;
  } else {
    return false;
  }
}

inline void fill_rppair(ReducedParticlePair &rpp, unsigned long id, double r21,
                        double q21) {
  rpp.id_ = id;
  rpp.r21_ = r21;
  rpp.q21_ = q21;
}

inline void fill_rppair(ReducedParticlePair &rpp, unsigned long id, double r21,
                        double q21, std::vector<long> repetitions) {
  rpp.id_ = id;
  rpp.r21_ = r21;
  rpp.q21_ = q21;
  rpp.repetitions_ = repetitions;
}

inline void fill_rppair(ReducedParticlePair &rpp, unsigned long id, double r21,
                        double q21, long irepeated) {
  rpp.id_ = id;
  rpp.r21_ = r21;
  rpp.q21_ = q21;
  rpp.repetitions_.push_back(irepeated);
}

struct reducedPairsComparison {
  bool operator()(const ReducedParticlePair &lx,
                  const ReducedParticlePair &rx) const {
    if (lx.r21_ != rx.r21_) {
      return (lx.r21_ < rx.r21_);
    }
    return (lx.q21_ < rx.q21_);
  }
};

inline void reducedparticlepair_tostream(const ReducedParticlePair &rpp,
                                         std::ostream &stream) {
  stream << rpp.id_ << '\t' << rpp.r21_ << '\t' << rpp.q21_ << '\t';
  // << "repeated = \t";
  for (auto ir : rpp.repetitions_) {
    stream << ir << ",";
  }
}

inline void reducedparticlepair_fromstream(ReducedParticlePair &rpp,
                                           std::istream &stream) {
  std::string sid_;
  std::string sr21_;
  std::string sq21_;

  std::string srep_;

  unsigned long id_;
  double r21_;  //!< Value for radius fraction r2/r1
  double q21_;  //!< Value for radius fraction r2/r1

  std::vector<long> repetitions_;

  stream >> sid_ >> sr21_ >> sq21_ >> srep_;

  id_ = boost::lexical_cast<unsigned long>(sid_);
  r21_ = std::stod(sr21_);
  q21_ = std::stod(sq21_);

  // split line
  std::vector<std::string> vssrep_;
  boost::split(vssrep_, srep_, boost::is_any_of(","));

  for (std::vector<std::string>::iterator it = vssrep_.begin();
       it < vssrep_.end(); ++it) {
    try {
      long trep = (boost::lexical_cast<long>(*it));
      repetitions_.push_back(trep);
    } catch (const boost::bad_lexical_cast &) {
      // std::cerr << "\n[ee] Error from lexical cast : " << *it;
      break;
    }
    // std::cout << "\t" << *it;
  }

  fill_rppair(rpp, id_, r21_, q21_, repetitions_);
}

class ContactPotential {
 public:
  ContactPotential(const unsigned long &id = 0, const double &potential = 0,
                   const short &n = 0)
      : id_(id), potential_(potential), n_(n) {}

  friend void fill_contact_potential(ContactPotential &cpot,
                                     const unsigned long &id,
                                     const double &potential, const short &n);

  void print(std::ostream &stream = std::cout) {
    contactpotential_tostream(*this, stream);
  }

  friend void contactpotential_tostream(const ContactPotential &cpot,
                                        std::ostream &stream);

  unsigned long id_;
  double potential_;
  short n_;
};

inline void fill_contact_potential(ContactPotential &cpot,
                                   const unsigned long &id,
                                   const double &potential, const short &n) {
  cpot.id_ = id;
  cpot.potential_ = potential;
  cpot.n_ = n;
}

inline void contactpotential_tostream(const ContactPotential &cpot,
                                      std::ostream &stream) {
  stream << '\n' << cpot.id_ << '\t' << cpot.potential_ << '\t' << cpot.n_;
}

struct find_id {
  unsigned long id;
  find_id(const unsigned long &id_) : id(id_) {}

  bool operator()(const ContactPotential &cpot) const {
    // return rpp.id_ = id
    // std::vector<long> rptd = rpp.repetitions_;
    // std::vector<long>::iterator it = std::find(rptd.begin(), rptd.end(), id);
    // return std::find(rptd.begin(), rptd.end(), id) != rptd.end();
    return (cpot.id_ == id);
  }
};

typedef std::vector<ContactPotential>::iterator ContactPotentialIterator;

class Enhancement {
 public:
  // Enhancement() {
  //}

  Enhancement(const darray &rarray_, const darray &qarray_,
              boost_array4d_ref efactor_, boost_array4d_ref cpotentials_,
              boost_array4d_ref bpotentials_, boost_array4d_ref rbarriers_,
              double eps_, src::severity_logger<severity_level> lg_)
      : rarray(rarray_),
        qarray(qarray_),
        efactor(efactor_),
        cpotentials(cpotentials_),
        bpotentials(bpotentials_),
        rbarriers(rbarriers_),
        eps(eps_),
        lg(lg_) {
    rarray_size = static_cast<unsigned short>(rarray.size());
    qarray_size = static_cast<unsigned short>(qarray.size());

    ncombs =
        ((rarray_size * qarray_size) * (rarray_size * qarray_size + 1)) / 2;

    particle_pairs.reserve(ncombs);

    potential_threshold = 0.0;

    nmin = 25;
    nmax = 2000;
    nstep = 5;

    AH = 20e-20;
  }

  Enhancement(const darray &rarray_, const darray &qarray_, double eps_,
              src::severity_logger<severity_level> lg_)
      : rarray(rarray_),
        qarray(qarray_),
        efactor(dummy),  // Hack, references must belong to an actual object
        cpotentials(dummy),
        bpotentials(dummy),
        rbarriers(dummy),
        eps(eps_),
        lg(lg_) {
    rarray_size = static_cast<unsigned short>(rarray.size());
    qarray_size = static_cast<unsigned short>(qarray.size());

    ncombs =
        ((rarray_size * qarray_size) * (rarray_size * qarray_size + 1)) / 2;

    particle_pairs.reserve(ncombs);

    potential_threshold = 0.0;

    nmin = 25;
    nmax = 2000;
    nstep = 5;
  }

  Enhancement(const darray &rarray_, const darray &qarray_,
              boost_array4d_ref efactor_, double eps_,
              src::severity_logger<severity_level> lg_)
      : rarray(rarray_),
        qarray(qarray_),
        efactor(efactor_),  // Hack, references must belong to an actual object
        cpotentials(dummy),
        bpotentials(dummy),
        rbarriers(dummy),
        eps(eps_),
        lg(lg_) {
    rarray_size = static_cast<unsigned short>(rarray.size());
    qarray_size = static_cast<unsigned short>(qarray.size());

    ncombs =
        ((rarray_size * qarray_size) * (rarray_size * qarray_size + 1)) / 2;

    particle_pairs.reserve(ncombs);

    potential_threshold = 0.0;

    nmin = 25;
    nmax = 2000;
    nstep = 5;
  }

  inline int write_particlepairs(std::string ppfilename) {
    std::ofstream ppstream(std::string(ppfilename + "_pp.dat"));

    ppstream << "#id\tl\tq\tm\tp\tr21\tq21\tnotswapd";
    // #pragma omp parallel for ordered//schedule(nonmonotonic:dynamic)
    for (unsigned int i = 0; i < particle_pairs.size(); ++i) {
      // #pragma omp ordered
      // {
      ppstream << '\n';
      particlepair_tostream(particle_pairs[i], ppstream);
      // }
    }
    ppstream.close();

    // Neutral particles
    std::ofstream npstream(std::string(ppfilename + "_np.dat"));

    npstream << "#id\tl\tq\tm\tp\tr21\tq21\tnotswapd";
    // #pragma omp parallel for ordered//schedule(nonmonotonic:dynamic)
    for (unsigned int i = 0; i < neutral_pairs.size(); ++i) {
      // #pragma omp ordered
      // {
      npstream << '\n';
      particlepair_tostream(neutral_pairs[i], npstream);
      // }
    }
    npstream.close();

    // Reduced pairs
    std::ofstream rpstream(std::string(ppfilename + "_rp.dat"));

    rpstream << "#id\tr21\tq21\trepetitions";
    // #pragma omp parallel for ordered//schedule(nonmonotonic:dynamic)
    for (unsigned int i = 0; i < reduced_pairs.size(); ++i) {
      // #pragma omp ordered
      // {
      rpstream << '\n';
      reducedparticlepair_tostream(reduced_pairs[i], rpstream);
      // }
    }
    rpstream.close();

    return 0;
  }

  inline int read_particlepairs(std::string ppfilename) {
    // Read particle pairs
    {
      // Clear vector
      particle_pairs.clear();

      std::ifstream ppstream(std::string(ppfilename + "_pp.dat"));

      std::string line;
      unsigned int nlines = 0;
      std::getline(ppstream, line);
      // if a comment first line isnt found, rewind the file
      if (line.find("#") == std::string::npos) {
        ppstream.clear();
        ppstream.seekg(0);
      }
      // read line by line
      while (std::getline(ppstream, line)) {
        ParticlePair pp;
        std::istringstream iss(line);
        particlepair_fromstream(pp, iss);
        particle_pairs.push_back(pp);
        ++nlines;
      }
      ppstream.close();
      BOOST_LOG_SEV(lg, info)
          << "Particle pair size = " << particle_pairs.size();
    }

    // Read particle pairs
    {
      // Clear vector
      neutral_pairs.clear();

      std::ifstream npstream(std::string(ppfilename + "_np.dat"));

      std::string line;
      unsigned int nlines = 0;
      std::getline(npstream, line);
      // if a comment first line isnt found, rewind the file
      if (line.find("#") == std::string::npos) {
        npstream.clear();
        npstream.seekg(0);
      }
      // read line by line
      while (std::getline(npstream, line)) {
        ParticlePair np;
        std::istringstream iss(line);
        particlepair_fromstream(np, iss);
        neutral_pairs.push_back(np);
        ++nlines;
      }
      npstream.close();
      BOOST_LOG_SEV(lg, info) << "Neutral pair size = " << neutral_pairs.size();
    }

    // Read reduced pairs
    {
      // Clear vector
      reduced_pairs.clear();

      std::ifstream rpstream(std::string(ppfilename + "_rp.dat"));

      std::string line;
      unsigned int nlines = 0;
      std::getline(rpstream, line);
      // if a comment first line isnt found, rewind the file
      if (line.find("#") == std::string::npos) {
        rpstream.clear();
        rpstream.seekg(0);
      }
      // read line by line
      while (std::getline(rpstream, line)) {
        ReducedParticlePair rp;
        std::istringstream iss(line);
        reducedparticlepair_fromstream(rp, iss);
        reduced_pairs.push_back(rp);
        ++nlines;
      }
      rpstream.close();
      BOOST_LOG_SEV(lg, info) << "Reduced pair size = " << reduced_pairs.size();
    }

    // for (auto rp: reduced_pairs){
    //   std::cout << std::endl;
    //   rp.print();
    // }
    if (reduced_pairs.size() > 0) {
      contact_potentials.resize(reduced_pairs.size());
    }

    return 0;
  }

  inline void compute_reducedpairs() {
    BOOST_LOG_SEV(lg, info) << "Computing particle pairs...";
    auto start = std::chrono::system_clock::now();
    //
    unsigned long idpp = 0;  // index for particle pairs
    unsigned long idnp = 0;  // index for neutral pairs
    // #pragma omp parallel for collapse(4) schedule(auto)
    for (unsigned int l = 0; l < rarray_size; ++l) {
      // iterate in charges particle 1
      for (unsigned int q = 0; q < qarray_size; ++q) {
        // iterate in radii particle 2
        for (unsigned int m = 0; m < rarray_size; ++m) {
          // iterate in charges particle 2
          for (unsigned int p = 0; p < qarray_size; ++p) {
            unsigned int ip1 = q * rarray_size + l;
            unsigned int ip2 = p * rarray_size + m;
            // avoid repetitions
            if ((p >= q) && (ip2 >= ip1)) {
              double r1 = rarray[l];
              double q1 = qarray[q];
              double r2 = rarray[m];
              double q2 = qarray[p];

              unsigned int lp1 = l;
              unsigned int qp1 = q;
              unsigned int mp2 = m;
              unsigned int pp2 = p;

              bool notswapd = true;

              // swap particles. We want to keep r21 < 1
              if (r2 > r1) {
                std::swap(r1, r2);
                std::swap(q1, q2);
                std::swap(lp1, mp2);
                std::swap(qp1, pp2);
                notswapd = false;
              }

              //double r21 = r2 / r1;
              double q21 = q2 / q1;

              // WARNING float comparison
              if (q1 == 0.0) {
                // q1=0, permute particle 1 with 2
                std::swap(r1, r2);
                std::swap(q1, q2);
                std::swap(lp1, mp2);
                std::swap(qp1, pp2);
                q21 = 0.0;
                notswapd = notswapd | false;
              }

              double r21 = r2 / r1;

              // if(q2==0.0){
              //   q21 = 0.0;
              // }

              // all neutrals case
              if ((q1 == 0.0) && (q2 == 0.0)) {
                // lp1 = l;  qp1 = q;
                // mp2 = m;  pp2 = p;
                // q21 = 0.0;
                // r21 = rarray[m]/rarray[l];
                ParticlePair neutralpair(idnp, lp1, qp1, mp2, pp2, r21, q21,
                                         notswapd);
                neutral_pairs.push_back(neutralpair);
                ++idnp;
                continue;  // break;
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

    std::chrono::duration<double> elapsed_seconds = end - start;
    BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();

    BOOST_LOG_SEV(lg, info) << "Pairs size : " << particle_pairs.size();
    BOOST_LOG_SEV(lg, info) << "Neutral pairs size : " << neutral_pairs.size();

    BOOST_LOG_SEV(lg, info) << "Total pairs size : "
                            << particle_pairs.size() + neutral_pairs.size();
    BOOST_LOG_SEV(lg, info) << "Total combinations : " << ncombs;

    BOOST_LOG_SEV(lg, info) << "Computing reduced pairs...";

    reduced_pairs.resize(particle_pairs.size());

    start = std::chrono::system_clock::now();
    for (unsigned int i = 0; i < particle_pairs.size(); ++i) {
      reduced_pairs[i].id_ = particle_pairs[i].id_;
      reduced_pairs[i].r21_ = particle_pairs[i].r21_;
      reduced_pairs[i].q21_ = particle_pairs[i].q21_;
    }

    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();

    start = std::chrono::system_clock::now();
    BOOST_LOG_SEV(lg, info) << "Sorting...";
    //__gnu_parallel::sort(reduced_pairs.begin(), reduced_pairs.end(),
    //reducedPairsComparison());
    std::sort(reduced_pairs.begin(), reduced_pairs.end(),
              reducedPairsComparison());
    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();

    // ******** new stuff
    // vector sorted with repeated
    std::vector<ReducedParticlePair> sorted_reducedpairs = reduced_pairs;
    // ******** end new stuff

    BOOST_LOG_SEV(lg, info) << "Erasing...";
    reduced_pairs.erase(std::unique(reduced_pairs.begin(), reduced_pairs.end()),
                        reduced_pairs.end());
    //__gnu_parallel::unique_copy(reduced_pairs.begin(), reduced_pairs.end(),
    //reduced_pairs);

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
    for (unsigned int irp = 0; irp < reduced_pairs.size(); ++irp) {
      for (long jpp = jstart; jpp < sorted_reducedpairs.size(); ++jpp) {
        if (reduced_pairs[irp] == sorted_reducedpairs[jpp]) {
          // #pragma omp critical
          // 	    {
          // one thread a at time
          reduced_pairs[irp].fill_repeated(sorted_reducedpairs[jpp].id_);
          // }
        } else {
          jstart = jpp;
          break;
        }
      }
    }

    /// end new stuff

    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();

    start = std::chrono::system_clock::now();
    // resize contact potentials
    if (reduced_pairs.size() > 0) {
      contact_potentials.resize(reduced_pairs.size());
    }

    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();

    BOOST_LOG_SEV(lg, info) << "Done computing pairs.";
  }

  inline 
  void compute_ipavdwpotential_contact(double cutoff, double vdw=1.0) {
    BOOST_LOG_SEV(lg, info) << "Computing IPA+vdW at contact...";
    auto start = std::chrono::system_clock::now();
    //
    unsigned long idpp = 0;  // index for particle pairs
    unsigned long idnp = 0;  // index for neutral pairs
    // #pragma omp parallel for collapse(4) schedule(auto)
    for (unsigned int l = 0; l < rarray_size; ++l) {
      // iterate in charges particle 1
      for (unsigned int q = 0; q < qarray_size; ++q) {
        // iterate in radii particle 2
        for (unsigned int m = 0; m < rarray_size; ++m) {
          // iterate in charges particle 2
          for (unsigned int p = 0; p < qarray_size; ++p) {
            double r1 = rarray[l];
            double q1 = qarray[q];
            double r2 = rarray[m];
            double q2 = qarray[p];
            double rt = r1+r2;
            eint::potential_ipavdw_funct pipavdw(r1, q1, r2, q2, eps, AH, cutoff, 1.0, vdw); 
            cpotentials[l][q][m][p] = pipavdw(rt);
          }
        }
      }
    }
    BOOST_LOG_SEV(lg, info) << "Done computing IPA+vdW potential contact.";
  }


  inline
  short particle_number(short l, short q, short L) {
    return l+L*q;
  }

  inline
  void mpcvdwpotential_barrier(short l, short q,
                               short m, short p,
                               double cutoff) {
    double r1 = rarray[l];
    double q1 = qarray[q];
    double r2 = rarray[m];
    double q2 = qarray[p];
    
    double rmin = r1+r2;
    double rmax = 100.0*rmin;
    
    double min = rmin;
    double max = rmax;

    short n = nterms4d[l][q][m][p];
    
    eint::force_mpcvdw_funct fmpcvdwfunct(r1, q1, r2, q2, eps, n, AH, cutoff, 1.0/*mpc*/, 1.0/*vdw*/);
    // force mpc at contact
    double fmpcvdw_rmin = fmpcvdwfunct(rmin);
    // Force at r max
    double fmpcvdw_rmax = fmpcvdwfunct(rmax);

    // checks if minimum exists
    if (fmpcvdw_rmin * fmpcvdw_rmax < 0.0) {
      // std::cerr << "\n[ii] Mixed phi_rt = " << fmpc_rmin << '\t' <<
      // fmpc_rmax;
      boost::uintmax_t bmax_iter = ROOT_MAXITER;
      tools::eps_tolerance<double> tol = ROOT_TOL;

      std::pair<double, double> pair_fmpcvdw;
      // try {
      pair_fmpcvdw = tools::toms748_solve(fmpcvdwfunct, min, max, fmpcvdw_rmin,
                                          fmpcvdw_rmax, tol, bmax_iter);
      if (bmax_iter > 990) {
        std::cerr << "\n ERROR max iter " << bmax_iter << "\n\n";
        std::terminate();
      }
      double rbarrier = 0.5 * (pair_fmpcvdw.first + pair_fmpcvdw.second);
      if (rbarrier >= min) {
        //*********************** USE SAME COEFFICIENTS OF FORCE
        //#pragma omp critical
        //{
          // mpc functor
          eint::potential_mpcvdw_funct pmpcvdwfunct(r1, q1, r2, q2, eps, n, AH, cutoff, 1.0/*mpc*/, 1.0/*vdw*/);

          // mpc at contact
          double pmpcvdw_rb = pmpcvdwfunct(rbarrier);

          bpotentials[l][q][m][p] = pmpcvdw_rb;
          rbarriers[l][q][m][p] = rbarrier;
          bpotentials[m][p][l][q] = bpotentials[l][q][m][p];
          rbarriers[m][p][l][q] = rbarriers[l][q][m][p];
        //}
      }
      else {
        std::cerr << "\n ERROR Negative rbarrier " << rbarrier << '\n';
        std::terminate();
      }
    }
  }                    
  inline
  void mpcvdwpotential_contact(short l, short q,
                               short m, short p,
                               double cutoff) {
    double r1 = rarray[l];
    double q1 = qarray[q];
    double r2 = rarray[m];
    double q2 = qarray[p];
    double rt = r1+r2;
    
    eint::potential_ipavdw_funct pipavdw(r1, q1, r2, q2, eps, AH,
                                          cutoff, 1.0, 1.0); 

    // ipa + vdw at contact
    double pipavdw_rt = pipavdw(rt);
    double pcomp = pipavdw_rt;

    unsigned short initer = 0;
    
    for (unsigned short n = nmin; n <= nmax; n += nstep, ++initer) {
      // mpc functor
      eint::potential_mpcvdw_funct pmpcvdwfunct(r1, q1, r2, q2, eps, n, 
                                                AH, cutoff, 1.0/*mpc*/, 1.0/*vdw*/);

      // mpc at contact
      double pmpcvdw_rt = pmpcvdwfunct(rt);

      double error_comp = max_pct_error(pcomp, pmpcvdw_rt);

      // std::cerr << '\n' << pcomp << '\t' << pmpcvdw_rt
      //                   << '\t' << cpotentials[l][q][m][p];

      if ((error_comp < MPC_ERROR) && (initer > 0)) {
        cpotentials[l][q][m][p] = pmpcvdw_rt;
        nterms4d[l][q][m][p] = n;
        cpotentials[m][p][l][q] = cpotentials[l][q][m][p];
        nterms4d[m][p][l][q] = nterms4d[l][q][m][p];
        // std::cerr << '\n' << pcomp << '\t' << pmpcvdw_rt
        //                 << '\t' << cpotentials[l][q][m][p]
        //                 << '\t' << n;
        // end iterations
        n = nmax;
      }
      else {
        if (n > nmax - nstep) {
          std::cerr << "\n[ww] Max iterations exceeded\n";
          cpotentials[l][q][m][p] = pmpcvdw_rt;
          nterms4d[l][q][m][p] = n;
          cpotentials[m][p][l][q] = cpotentials[l][q][m][p];
          nterms4d[m][p][l][q] = nterms4d[l][q][m][p];
        }
      }
      pcomp = pmpcvdw_rt;
    }
  }

  template <typename PotentialFunctor>
  void iterate_symmetric(short L, short Q, PotentialFunctor func,
                         double cutoff) {
    // n differents particles with r and q
    short n = L*Q;

    // vector of particle numbers -> pair l, q
    vp particle(n);
    for(short l=0; l<L; ++l) {
      for(short q=0; q<Q; ++q) {
          pp p = pp(l,q);
          short pnumber = particle_number(l, q, L);
          particle[pnumber] = p;
      }
    }

    // number of unique potentials to compute - non repeated combinations
    // of particles    
    unsigned comb = n*(n+1)/2;

    // vector to store the combinations
    vp vps(comb);
    
    unsigned count = 0;
    for(short i=0; i<n; ++i) {
      for(short k=i; k<n; ++k) {
          pp p = pp(i, k);
          // store the combinations particle number i , particle number j
          vps[count] = p;
          count++;
      }
    }
    BOOST_LOG_SEV(lg, info) << "\n total length " << vps.size();
    count = 0;

    #pragma omp parallel for
    for(unsigned k=0; k<vps.size(); k++) {

      short i = vps[k].first;// particle i
      short j = vps[k].second;// particle j

      count++;

      // get indices (l,q), (m,r) for particles i and j
      short l = particle[i].first;
      short q = particle[i].second;
      short m = particle[j].first;
      short r = particle[j].second;

      // apply function to indices and cutoff
      (this->*func)(l, q, m, r, cutoff);
    }
  }


  inline
  void compute_mpcvdwpotential_contact_sym(double cutoff) {
    BOOST_LOG_SEV(lg, info) << "Computing MPC+vdW at contact symmetric...";

    iterate_symmetric(rarray_size, qarray_size, &Enhancement::mpcvdwpotential_contact, cutoff);

  }

  inline
  void compute_mpcvdwpotential_barrier_sym(double cutoff) {
    BOOST_LOG_SEV(lg, info) << "Computing MPC+vdW barrier symmetric...";
    // pass the reference of the method to iterate_symmetric
    iterate_symmetric(rarray_size, qarray_size, &Enhancement::mpcvdwpotential_barrier, cutoff);

  }

  inline
  void compute_mpcvdwpotential_contact_bf(double cutoff) {

    BOOST_LOG_SEV(lg, info) << "Computing MPC+vdW at contact...";

    #pragma omp parallel for ordered collapse(4) schedule(auto)
    for (unsigned int l = 0; l < rarray_size; ++l) {
      // iterate in charges particle 1
      for (unsigned int q = 0; q < qarray_size; ++q) {
        // iterate in radii particle 2
        for (unsigned int m = 0; m < rarray_size; ++m) {
          // iterate in charges particle 2
          for (unsigned int p = 0; p < qarray_size; ++p) {
            double r1 = rarray[l];
            double q1 = qarray[q];
            double r2 = rarray[m];
            double q2 = qarray[p];
            double rt = r1+r2;
            
            eint::potential_ipavdw_funct pipavdw(r1, q1, r2, q2, eps, AH,
                                                 cutoff, 1.0, 1.0); 

            // ipa + vdw at contact
            double pipavdw_rt = pipavdw(rt);
            double pcomp = pipavdw_rt;

            unsigned short initer = 0;
            
            for (unsigned short n = nmin; n <= nmax; n += nstep, ++initer) {
              // mpc functor
              eint::potential_mpcvdw_funct pmpcvdwfunct(r1, q1, r2, q2, eps, n, 
                                                        AH, cutoff, 1.0/*mpc*/, 1.0/*vdw*/);

              // mpc at contact
              double pmpcvdw_rt = pmpcvdwfunct(rt);

              double error_comp = max_pct_error(pcomp, pmpcvdw_rt);

              // std::cerr << '\n' << pcomp << '\t' << pmpcvdw_rt
              //                   << '\t' << cpotentials[l][q][m][p];

              if ((error_comp < MPC_ERROR) && (initer > 0)) {
                cpotentials[l][q][m][p] = pmpcvdw_rt;
                nterms4d[l][q][m][p] = n;
                // std::cerr << '\n' << pcomp << '\t' << pmpcvdw_rt
                //                 << '\t' << cpotentials[l][q][m][p]
                //                 << '\t' << n;
                // end iterations
                n = nmax;
              }
              else {
                if (n > nmax - nstep) {
                  std::cerr << "\n[ww] Max iterations exceeded\n";
                  cpotentials[l][q][m][p] = pmpcvdw_rt;
                  nterms4d[l][q][m][p] = n;
                }
              }
              pcomp = pmpcvdw_rt;
            }            
          }
        }
      }
    }
    BOOST_LOG_SEV(lg, info) << "Done computing MPC+vdW potential contact.";
  }

  // resize nterms4d grid for mpc
  inline
  void set_nterms(bgrid4d grid4) {
    nterms4d.resize(grid4);
  }

  inline
  void compute_mpcvdwpotential_barrier_bf(double cutoff) {
    BOOST_LOG_SEV(lg, info) << "Computing MPC+vdW barriers...";
    auto start = std::chrono::system_clock::now();
    //
    #pragma omp parallel for ordered collapse(4) schedule(auto)
    for (unsigned int l = 0; l < rarray_size; ++l) {
      // iterate in charges particle 1
      for (unsigned int q = 0; q < qarray_size; ++q) {
        // iterate in radii particle 2
        for (unsigned int m = 0; m < rarray_size; ++m) {
          // iterate in charges particle 2
          for (unsigned int p = 0; p < qarray_size; ++p) {
            double r1 = rarray[l];
            double q1 = qarray[q];
            double r2 = rarray[m];
            double q2 = qarray[p];
            
            double rmin = r1+r2;
            double rmax = 100.0*rmin;
            
            double min = rmin;
            double max = rmax;

            short n = nterms4d[l][q][m][p];
            
            eint::force_mpcvdw_funct fmpcvdwfunct(r1, q1, r2, q2, eps, n, AH, cutoff, 1.0/*mpc*/, 1.0/*vdw*/);
            // force mpc at contact
            double fmpcvdw_rmin = fmpcvdwfunct(rmin);
            // Force at r max
            double fmpcvdw_rmax = fmpcvdwfunct(rmax);

            // checks if minimum exists
            if (fmpcvdw_rmin * fmpcvdw_rmax < 0.0) {
              // std::cerr << "\n[ii] Mixed phi_rt = " << fmpc_rmin << '\t' <<
              // fmpc_rmax;
              boost::uintmax_t bmax_iter = ROOT_MAXITER;
              tools::eps_tolerance<double> tol = ROOT_TOL;

              std::pair<double, double> pair_fmpcvdw;
              // try {
              pair_fmpcvdw = tools::toms748_solve(fmpcvdwfunct, min, max, fmpcvdw_rmin,
                                                  fmpcvdw_rmax, tol, bmax_iter);
              if (bmax_iter > 990) {
                std::cerr << "\n ERROR max iter " << bmax_iter << "\n\n";
                std::terminate();
              }
              double rbarrier = 0.5 * (pair_fmpcvdw.first + pair_fmpcvdw.second);
              if (rbarrier >= min) {
                //*********************** USE SAME COEFFICIENTS OF FORCE
                //#pragma omp critical
                //{
                  // mpc functor
                  eint::potential_mpcvdw_funct pmpcvdwfunct(r1, q1, r2, q2, eps, n, AH, cutoff, 1.0/*mpc*/, 1.0/*vdw*/);

                  // mpc at contact
                  double pmpcvdw_rb = pmpcvdwfunct(rbarrier);

                  bpotentials[l][q][m][p] = pmpcvdw_rb;
                  rbarriers[l][q][m][p] = rbarrier;
                //}
              }
              else {
                std::cerr << "\n ERROR Negative rbarrier " << rbarrier << '\n';
                std::terminate();
              }
            }
          }
        }
      }
    }

    BOOST_LOG_SEV(lg, info) << "Done computing MPC+vdW potential barrier.";
  }

  inline
  void compute_ipavdwpotential_barrier(double cutoff, double vdw=1.0) {
    BOOST_LOG_SEV(lg, info) << "Computing IPA+vdW barrier...";
    auto start = std::chrono::system_clock::now();
    //
    // #pragma omp parallel for collapse(4) schedule(auto)
    for (unsigned int l = 0; l < rarray_size; ++l) {
      // iterate in charges particle 1
      for (unsigned int q = 0; q < qarray_size; ++q) {
        // iterate in radii particle 2
        for (unsigned int m = 0; m < rarray_size; ++m) {
          // iterate in charges particle 2
          for (unsigned int p = 0; p < qarray_size; ++p) {
            double r1 = rarray[l];
            double q1 = qarray[q];
            double r2 = rarray[m];
            double q2 = qarray[p];
            
            double rmin = r1+r2;
            double rmax = 100.0*rmin;
            
            double min = rmin;
            double max = rmax;            

            eint::force_ipavdw_funct fipavdw(r1, q1, r2, q2, eps, AH, cutoff, 1.0, vdw);
            // force ipa at contact
            double fipa_rmin = fipavdw(rmin);
            // Force at r max
            double fipa_rmax = fipavdw(rmax);

            // checks if minimum exists
            if (fipa_rmin * fipa_rmax < 0.0) {
              // std::cerr << "\n[ii] Mixed phi_rt = " << fipa_rmin << '\t' <<
              // fipa_rmax;
              boost::uintmax_t bmax_iter = ROOT_MAXITER;
              tools::eps_tolerance<double> tol = ROOT_TOL;

              std::pair<double, double> pair_fipa;
              //try {
              pair_fipa = tools::toms748_solve(fipavdw, min, max, fipa_rmin,
                                                 fipa_rmax, tol, bmax_iter);
              // } catch (const std::exception &exc) {
              //   std::cerr << '\n' << exc.what() << '\n';
              //   // std::terminate();
              // std::cerr << "\n[ee] No barrier\n";
              // }
              if (bmax_iter > 990) {
                std::cerr << "\n ERROR max iter " << bmax_iter << "\n\n";
                std::terminate();
              }
              double rbarrier = 0.5 * (pair_fipa.first + pair_fipa.second);
              if (rbarrier >= min) {
                //*********************** USE SAME COEFFICIENTS OF FORCE
                // ipa functor
                eint::potential_ipavdw_funct pipavdw(r1, q1, r2, q2, eps, AH, cutoff, 1.0, 1.0); 
                // ipa at contact
                double pipa_rbarrier = pipavdw(rbarrier);
                // the barrier has the index of the id of reduced pairs
      

                  bpotentials[l][q][m][p] = pipa_rbarrier;
                  rbarriers[l][q][m][p] = rbarrier;
              }
              else {
                std::cerr << "\n ERROR Negative rbarrier " << rbarrier << '\n';
                std::terminate();
              }
            }
          }
        }
      }
    }

    BOOST_LOG_SEV(lg, info) << "Done computing IPA+vdW potential barrier.";
  }

  // compute bruteforce
  // iterates in the whole grid
  inline
  int compute_bruteforce() {
    auto start = std::chrono::system_clock::now();
    #pragma omp parallel for collapse(4) ordered
    for (unsigned int l = 0; l < rarray_size; ++l) {
      // iterate in charges particle 1
      for (unsigned int q = 0; q < qarray_size; ++q) {
        // iterate in radii particle 2
        for (unsigned int m = 0; m < rarray_size; ++m) {
          // iterate in charges particle 2
          for (unsigned int p = 0; p < qarray_size; ++p) {
            #pragma omp critical 
            {
            double phimin = cpotentials[l][q][m][p];
            double phimax = bpotentials[l][q][m][p];

            // checks if a barrier exists over potential_threshold
            // FIXME change to 0
            // if (fabs(phimax) > potential_threshold) {
            if (fabs(phimax) != 0.0) {
              // double phimax = bpotentials[l][q][m][p];
              double eta = eta_barrier(phimin, phimax);

              // In the hybrid approach eta can be negative
              // in the case of phimin > phimax
              
              if (eta < 0.0) {
                if (fabs(eta) > ETA_TOL) {
                  BOOST_LOG_SEV(lg, warning) << "Negative enhancement factor : "
                                           << eta;
                  BOOST_LOG_SEV(lg, warning) << "phimin : " << phimin;
                  BOOST_LOG_SEV(lg, warning) << "phimax : " << phimax;
                  //BOOST_LOG_SEV(lg, warning) << "Abort.";
                  //std::terminate();
                  if (phimin > 0.0) {
                    eta = eta_repulsive(phimin);
                    BOOST_LOG_SEV(lg, warning) << "choosing repulsive phimin, eta: "
                                             << eta;
                  }
                  else {
                    eta = eta_attractive(phimin);
                    BOOST_LOG_SEV(lg, warning) << "choosing attractive phimax, eta: "
                                             << eta;
                  }
                }
                else {
                  BOOST_LOG_SEV(lg, warning) << "eta is negative but less than tolerance";
                  eta = 0.0;
                }
              }
              //#pragma omp atomic write
              efactor[l][q][m][p] = eta;
            } 
            else {
              if (phimin > 0.0) {
                double eta = eta_repulsive(phimin);
                //#pragma omp atomic write
                if (fabs(eta) < ETA_TOL) {
                  BOOST_LOG_SEV(lg, warning) << "eta is less than tolerance";
                  eta = 0.0;
                }
                efactor[l][q][m][p] = eta;
              }
              else {
                double eta = eta_attractive(phimin);
                //#pragma omp atomic write
                if (fabs(eta) < ETA_TOL) {
                  BOOST_LOG_SEV(lg, warning) << "eta is less than tolerance";
                  eta = 0.0;
                }
                efactor[l][q][m][p] = eta;
              }
            }
            }
          }
        }
      }
    }   
    //
    auto end = std::chrono::system_clock::now();
    auto elapsed_seconds = end-start;
    BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();
    return 0;
  }

  inline void compute_coulombpotential_contact() {
    BOOST_LOG_SEV(lg, info)
        << "Computing Coulomb potential at contact for n pairs : "
        << reduced_pairs.size();
#pragma omp parallel for ordered  // schedule(nonmonotonic:dynamic)
    for (unsigned int i = 0; i < reduced_pairs.size(); ++i) {
      unsigned long id = reduced_pairs[i].id_;
      double r21 = reduced_pairs[i].r21_;
      double q21 = reduced_pairs[i].q21_;

      double rt = 1.0 + r21;

      eint::potential_coulomb_funct pcoulfunct(r21, q21);

      // coulomb at contact
      double pcoul_rt = pcoulfunct(rt);

#pragma omp critical
      { fill_contact_potential(contact_potentials[i], id, pcoul_rt, 0); }
    }

    barrier_potentials.resize(0);
    rbarrier_array.resize(0);

    BOOST_LOG_SEV(lg, info) << "Done computing Coulomb potential at contact.";
  }

  inline void compute_coulombvdwpotential_contact(double cutoff) {
    BOOST_LOG_SEV(lg, info)
        << "Computing Coulomb+vdW potential at contact for n pairs : ";
    auto start = std::chrono::system_clock::now();
    //
    unsigned long idpp = 0;  // index for particle pairs
    unsigned long idnp = 0;  // index for neutral pairs
    // #pragma omp parallel for collapse(4) schedule(auto)
    for (unsigned int l = 0; l < rarray_size; ++l) {
      // iterate in charges particle 1
      for (unsigned int q = 0; q < qarray_size; ++q) {
        // iterate in radii particle 2
        for (unsigned int m = 0; m < rarray_size; ++m) {
          // iterate in charges particle 2
          for (unsigned int p = 0; p < qarray_size; ++p) {
            double r1 = rarray[l];
            double q1 = qarray[q];
            double r2 = rarray[m];
            double q2 = qarray[p];
            double rt = r1+r2;
            eint::potential_coulombvdw_funct pcoulvdw(r1, q1, r2, q2, eps, AH, cutoff, 1.0, 1.0); 
            cpotentials[l][q][m][p] = pcoulvdw(rt);
          }
        }
      }
    }
    BOOST_LOG_SEV(lg, info) << "Done computing IPA+vdW potential contact.";
  }

  inline
  void compute_coulombvdwpotential_barrier(double cutoff) {
    BOOST_LOG_SEV(lg, info) << "Computing Coulomb+vdW barrier...";
    auto start = std::chrono::system_clock::now();
    //
    // #pragma omp parallel for collapse(4) schedule(auto)
    for (unsigned int l = 0; l < rarray_size; ++l) {
      // iterate in charges particle 1
      for (unsigned int q = 0; q < qarray_size; ++q) {
        // iterate in radii particle 2
        for (unsigned int m = 0; m < rarray_size; ++m) {
          // iterate in charges particle 2
          for (unsigned int p = 0; p < qarray_size; ++p) {
            double r1 = rarray[l];
            double q1 = qarray[q];
            double r2 = rarray[m];
            double q2 = qarray[p];
            
            double rmin = r1+r2;
            double rmax = 100.0*rmin;
            
            double min = rmin;
            double max = rmax;            

            eint::force_coulombvdw_funct fcoulombvdw(r1, q1, r2, q2, eps, AH, cutoff, 1.0, 1.0);
            // force coulomb at contact
            double fcoulomb_rmin = fcoulombvdw(rmin);
            // Force at r max
            double fcoulomb_rmax = fcoulombvdw(rmax);

            // checks if minimum exists
            if (fcoulomb_rmin * fcoulomb_rmax < 0.0) {
              boost::uintmax_t bmax_iter = ROOT_MAXITER;
              tools::eps_tolerance<double> tol = ROOT_TOL;

              std::pair<double, double> pair_fcoulomb;
              
              pair_fcoulomb = tools::toms748_solve(fcoulombvdw, min, max, fcoulomb_rmin,
                                                 fcoulomb_rmax, tol, bmax_iter);
              
              if (bmax_iter > 990) {
                std::cerr << "\n ERROR max iter " << bmax_iter << "\n\n";
                std::terminate();
              }
              double rbarrier = 0.5 * (pair_fcoulomb.first + pair_fcoulomb.second);
              if (rbarrier >= min) {
                //*********************** USE SAME COEFFICIENTS OF FORCE
                // coulomb functor
                eint::potential_coulombvdw_funct pcoulombvdw(r1, q1, r2, q2, eps, AH, cutoff, 1.0, 1.0); 
                // coulomb at contact
                double pcoulomb_rbarrier = pcoulombvdw(rbarrier);
                // the barrier has the index of the id of reduced pairs
      

                  bpotentials[l][q][m][p] = pcoulomb_rbarrier;
                  rbarriers[l][q][m][p] = rbarrier;
              }
              else {
                std::cerr << "\n ERROR Negative rbarrier " << rbarrier << '\n';
                std::terminate();
              }
            }
          }
        }
      }
    }

    BOOST_LOG_SEV(lg, info) << "Done computing IPA+vdW potential barrier.";
  }

  inline void compute_ipapotential_contact() {
    BOOST_LOG_SEV(lg, info)
        << "Computing ipa potential at contact for n pairs : "
        << reduced_pairs.size();
#pragma omp parallel for ordered  // schedule(nonmonotonic:dynamic)
    for (unsigned int i = 0; i < reduced_pairs.size(); ++i) {
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
    }

    barrier_potentials.reserve(contact_potentials.size());
    rbarrier_array.reserve(contact_potentials.size());

    BOOST_LOG_SEV(lg, info) << "Done computing ipa potential at contact.";
  }

  inline void compute_ipapotential_barrier() {
    BOOST_LOG_SEV(lg, info) << "Computing ipa potential barrier for n pairs : "
                            << reduced_pairs.size();
#pragma omp \
    parallel for  // shared(reduced_pairs) schedule(nonmonotonic:dynamic)
    for (unsigned int i = 0; i < contact_potentials.size(); ++i) {
      unsigned int index = i;  // contact_potentials[i].id_;

      double r21 = reduced_pairs[index].r21_;
      double q21 = reduced_pairs[index].q21_;

      double rmin = 1.0 + r21;
      double rmax = 100.0 * rmin;

      double min = rmin;
      double max = rmax;

      eint::force_ipa_funct fipafunct(r21, q21, eps);
      // force ipa at contact
      double fipa_rmin = fipafunct(rmin);
      // Force at r max
      double fipa_rmax = fipafunct(rmax);

      // checks if minimum exists
      //	#pragma omp ordered
      if (fipa_rmin * fipa_rmax < 0.0) {
        // std::cerr << "\n[ii] Mixed phi_rt = " << fipa_rmin << '\t' <<
        // fipa_rmax;
        boost::uintmax_t bmax_iter = ROOT_MAXITER;
        tools::eps_tolerance<double> tol = ROOT_TOL;

        std::pair<double, double> pair_fipa;
        try {
          pair_fipa = tools::toms748_solve(fipafunct, min, max, fipa_rmin,
                                           fipa_rmax, tol, bmax_iter);
        } catch (const std::exception &exc) {
          std::cerr << '\n' << exc.what() << '\n';
          // std::terminate();
          std::cerr << "\n[ee] No barrier\n";
        }
        if (bmax_iter > 990) {
          std::cerr << "\n ERROR max iter " << bmax_iter << "\n\n";
          std::terminate();
        }
        double rbarrier = 0.5 * (pair_fipa.first + pair_fipa.second);
        if (rbarrier >= min) {
          //*********************** USE SAME COEFFICIENTS OF FORCE
          // ipa functor
          eint::potential_ipa_funct pipafunct(r21, q21, eps);
          // ipa at contact
          double pipa_rbarrier = pipafunct(rbarrier);
          // the barrier has the index of the id of reduced pairs
#pragma omp critical
          {
            ContactPotential barrierpot(index, pipa_rbarrier);
            barrier_potentials.push_back(barrierpot);
            rbarrier_array.push_back(rbarrier);
          }
        } else {
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

    BOOST_LOG_SEV(lg, info)
        << "Barrier potentials size : " << barrier_potentials.size();
    BOOST_LOG_SEV(lg, info) << "Done computing ipa potential barrier.";
  }

  inline void compute_mpcpotential_contact() {
    BOOST_LOG_SEV(lg, info)
        << "Computing mpc potential at contact for n pairs : "
        << reduced_pairs.size();
    // WARNING FIXME
    // std::ofstream outfile("pot1.dat");
    // outfile << "#n\tr21\tq21\tmpc\tpipa\terror_pot\terror\tmsg\n";

#pragma omp parallel for  // schedule(nonmonotonic:dynamic)
    for (unsigned int i = 0; i < reduced_pairs.size(); ++i) {
      unsigned long id = reduced_pairs[i].id_;
      double r21 = reduced_pairs[i].r21_;
      double q21 = reduced_pairs[i].q21_;

      double rt = 1.0 + r21;

      eint::potential_ipa_funct pipafunct(r21, q21, eps);

      // // ipa at contact
      double pipa_rt = pipafunct(rt);

      double pcomp = pipa_rt;
      unsigned int initer = 0;

      for (short n = nmin; n <= nmax; n += nstep, ++initer) {
        // mpc functor
        eint::potential_mpc_funct pmpcfunct(r21, q21, eps, n);

        // mpc at contact
        double pmpc_rt = pmpcfunct(rt);

        double error_comp = max_pct_error(pcomp, pmpc_rt);
        // double error_pipa = max_pct_error(pipa_rt, pmpc_rt);

        // cerr << "\n[ii] n = " << n << "\t err = " << error_comp;
        if ((error_comp < MPC_ERROR) && (initer > 0)) {
          // #pragma omp critical
          // 	    {
          fill_contact_potential(contact_potentials[i], id, pmpc_rt, n);
          // outfile << n
          // 	      << '\t' << r21 << '\t' << q21
          // 	      << '\t' << pmpc_rt
          // 	      << '\t' << pipa_rt
          // 	      << '\t' << error_pipa
          // 	      << '\t' << error_comp << '\n';

          n = nmax;

        } else {
          if (n > nmax - nstep) {
            std::cerr << "\n[ww] Max iterations exceeded\n";
            fill_contact_potential(contact_potentials[i], i, pmpc_rt, n);
          }
        }
        pcomp = pmpc_rt;
      }
    }

    //      outfile.close();

    barrier_potentials.reserve(contact_potentials.size());
    rbarrier_array.reserve(contact_potentials.size());

    BOOST_LOG_SEV(lg, info)
        << "Error potentials size : " << error_potentials.size();
    BOOST_LOG_SEV(lg, info) << "Done computing MPC potential at contact.";
  }

  inline double eta_attractive(double phimin) {
    double kt = Kboltz * temperature;
    return 1.0 - phimin / kt;
  }

  inline double eta_repulsive(double phimin) {
    double kt = Kboltz * temperature;
    return exp(-phimin / kt);
  }

  inline double eta_barrier(double phimin, double phimax) {
    double kt = Kboltz * temperature;
    return exp(-phimax / kt) * (1.0 + (phimax - phimin) / kt);
  }

  inline void compute_mpcpotential_barrier() {
    BOOST_LOG_SEV(lg, info) << "Computing mpc potential barrier for n pairs : "
                            << contact_potentials.size();
#pragma omp \
    parallel for  // shared(reduced_pairs) schedule(nonmonotonic:dynamic)
    for (unsigned int i = 0; i < contact_potentials.size(); ++i) {
      unsigned int index = i;  // contact_potentials[i].id_;
      unsigned int nterms = contact_potentials[i].n_;

      double r21 = reduced_pairs[index].r21_;
      double q21 = reduced_pairs[index].q21_;

      double rmin = 1.0 + r21;
      double rmax = 100.0 * rmin;

      double min = rmin;
      double max = rmax;

      eint::force_mpc_funct fmpcfunct(r21, q21, eps, nterms);
      // force mpc at contact
      double fmpc_rmin = fmpcfunct(rmin);
      // Force at r max
      double fmpc_rmax = fmpcfunct(rmax);

      // checks if minimum exists
      if (fmpc_rmin * fmpc_rmax < 0.0) {
        // std::cerr << "\n[ii] Mixed phi_rt = " << fmpc_rmin << '\t' <<
        // fmpc_rmax;
        boost::uintmax_t bmax_iter = ROOT_MAXITER;
        tools::eps_tolerance<double> tol = ROOT_TOL;

        std::pair<double, double> pair_fmpc;
        // try {
        pair_fmpc = tools::toms748_solve(fmpcfunct, min, max, fmpc_rmin,
                                         fmpc_rmax, tol, bmax_iter);
        // pair_fmpc = tools::bisect(fmpcfunct, min, max, tol, bmax_iter);
        // }
        // catch(const std::exception& exc) {
        //   std::cerr << '\n' << exc.what() << '\n';
        //   std::terminate();
        // }
        if (bmax_iter > 990) {
          std::cerr << "\n ERROR max iter " << bmax_iter << "\n\n";
          std::terminate();
        }
        double rbarrier = 0.5 * (pair_fmpc.first + pair_fmpc.second);
        if (rbarrier >= min) {
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
        } else {
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

    BOOST_LOG_SEV(lg, info)
        << "Barrier potentials size : " << barrier_potentials.size();
    BOOST_LOG_SEV(lg, info) << "Done computing mpc potential barrier.";
  }

  inline
  void expand_barriers() {
    auto start = std::chrono::system_clock::now();

    BOOST_LOG_SEV(lg, info) << "Size for array efactor : " << efactor.size();

    BOOST_LOG_SEV(lg, info) << "Enhancement: expanding pair potentials";

// iterate in reduced_pairs
#pragma omp parallel for
    for (unsigned int irp = 0; irp < reduced_pairs.size(); ++irp) {
      // this is the index of reduced pairs, i.e., non repeated pair
      unsigned int index = reduced_pairs[irp].id_;

      // contact potential for this reduced pairhas the same ordering
      // of reduced pairs
      double contact_potential = contact_potentials[irp].potential_;

      // get repeated combinations (of particle_pairs)
      std::vector<long> reps = reduced_pairs[irp].repetitions_;

      // find if index is in barrier_potentials, i.e., if a potential
      // barrier exists
      ContactPotentialIterator barrier_it = std::find_if(
          barrier_potentials.begin(), barrier_potentials.end(), find_id(irp));

      // iterate in vector of repeated combinations
      for (unsigned int jrep = 0; jrep < reps.size(); ++jrep) {
        long rep_index = reps[jrep];
        short l = particle_pairs[rep_index].l_;
        short q = particle_pairs[rep_index].q_;
        short m = particle_pairs[rep_index].m_;
        short p = particle_pairs[rep_index].p_;

        double tr1 = rarray[l];  //(!notswapd ? rarray[l] : rarray[m]);
        // std::cerr << "\nn q is " << q;
        double tq1 = qarray[q];  //(!notswapd ? qarray[q] : qarray[p]);
        double potprefactor = tq1 * tq1 / tr1;
        double phimin = potprefactor * contact_potential;
        // std::cerr << "\n\n";

        // WARNING HACK do nothing for contact potentials
        //cpotentials[m][p][l][q] = phimin;
        //cpotentials[l][q][m][p] = phimin;
        // we have potential barrier
        // ** if barrier compute eta
        if (barrier_it != barrier_potentials.end()) {
          double phimax = potprefactor * (*barrier_it).potential_;
          // update potentials
          bpotentials[m][p][l][q] = phimax;
          bpotentials[l][q][m][p] = phimax;
          unsigned int idr = barrier_it - barrier_potentials.begin();
          rbarriers[m][p][l][q] = rbarrier_array[idr] * tr1;
          rbarriers[l][q][m][p] = rbarrier_array[idr] * tr1;
        }

        // std::cout << std::endl << l << '\t' << q << '\t' << m << '\t' << p <<
        // '\t' << efactor[l][q][m][p];
      }
    }

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    BOOST_LOG_SEV(lg, info) << "Done expanding pair potentials";
  }

  inline void compute_enhancement_factor() {
    // for (unsigned int l=0; l<rarray_size; ++l) {
    // 	// iterate in charges particle 1
    // 	for (unsigned int q=0; q<qarray_size; ++q) {
    // 	  // iterate in radii particle 2
    // 	  for (unsigned int m=0; m<rarray_size; ++m) {
    // 	    // iterate in charges particle 2
    // 	    for (unsigned int p=0; p<qarray_size; ++p) {
    // 	      efactor[l][q][m][p] = 0.0;
    // 	      cpotentials[l][q][m][p] = 0.0;
    // 	      bpotentials[l][q][m][p] = 0.0;
    // 	    }
    // 	  }
    // 	}
    // }

    auto start = std::chrono::system_clock::now();

    BOOST_LOG_SEV(lg, info) << "Size for array efactor : " << efactor.size();

    BOOST_LOG_SEV(lg, info) << "Enhancement: expanding pair potentials";

// iterate in reduced_pairs
#pragma omp parallel for
    for (unsigned int irp = 0; irp < reduced_pairs.size(); ++irp) {
      // this is the index of reduced pairs, i.e., non repeated pair
      unsigned int index = reduced_pairs[irp].id_;

      // contact potential for this reduced pairhas the same ordering
      // of reduced pairs
      double contact_potential = contact_potentials[irp].potential_;

      // get repeated combinations (of particle_pairs)
      std::vector<long> reps = reduced_pairs[irp].repetitions_;

      // find if index is in barrier_potentials, i.e., if a potential
      // barrier exists
      ContactPotentialIterator barrier_it = std::find_if(
          barrier_potentials.begin(), barrier_potentials.end(), find_id(irp));

      // iterate in vector of repeated combinations
      for (unsigned int jrep = 0; jrep < reps.size(); ++jrep) {
        long rep_index = reps[jrep];
        short l = particle_pairs[rep_index].l_;
        short q = particle_pairs[rep_index].q_;
        short m = particle_pairs[rep_index].m_;
        short p = particle_pairs[rep_index].p_;

        double tr1 = rarray[l];  //(!notswapd ? rarray[l] : rarray[m]);
        // std::cerr << "\nn q is " << q;
        double tq1 = qarray[q];  //(!notswapd ? qarray[q] : qarray[p]);
        double potprefactor = tq1 * tq1 / tr1;
        double phimin = potprefactor * contact_potential;
        // std::cerr << "\n\n";

        cpotentials[m][p][l][q] = phimin;
        cpotentials[l][q][m][p] = phimin;
        // we have potential barrier
        // ** if barrier compute eta
        if (barrier_it != barrier_potentials.end()) {
          double phimax = potprefactor * (*barrier_it).potential_;
          // update potentials
          bpotentials[m][p][l][q] = phimax;
          bpotentials[l][q][m][p] = phimax;
          unsigned int idr = barrier_it - barrier_potentials.begin();
          rbarriers[m][p][l][q] = rbarrier_array[idr] * tr1;
          rbarriers[l][q][m][p] = rbarrier_array[idr] * tr1;
        }

        // std::cout << std::endl << l << '\t' << q << '\t' << m << '\t' << p <<
        // '\t' << efactor[l][q][m][p];
      }
    }

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    BOOST_LOG_SEV(lg, info) << "Done expanding pair potentials";
    BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();

    BOOST_LOG_SEV(lg, info) << "Computing Enhancement Factor";
    start = std::chrono::system_clock::now();

#pragma omp parallel for
    for (unsigned int l = 0; l < rarray_size; ++l) {
      // iterate in charges particle 1
      for (unsigned int q = 0; q < qarray_size; ++q) {
        // iterate in radii particle 2
        for (unsigned int m = 0; m < rarray_size; ++m) {
          // iterate in charges particle 2
          for (unsigned int p = 0; p < qarray_size; ++p) {
            double phimin = cpotentials[l][q][m][p];
            double phimax = bpotentials[l][q][m][p];

            // checks if a barrier exists over potential_threshold
            // FIXME change to 0
            // if (fabs(phimax) > potential_threshold) {
            if (fabs(phimax) != 0.0) {
              // double phimax = bpotentials[l][q][m][p];
              double eta = eta_barrier(phimin, phimax);

              // In the hybrid approach eta can be negative
              // in the case of phimin > phimax
              if (eta < 0.0) {
                BOOST_LOG_SEV(lg, warning) << "Negative enhancement factor";
                BOOST_LOG_SEV(lg, warning) << "phimin : " << phimin;
                BOOST_LOG_SEV(lg, warning) << "phimax : " << phimax;
                if (phimin > 0.0) {
                  eta = eta_repulsive(phimin);
                } else {
                  eta = eta_attractive(phimin);
                }
              }
              //#pragma omp atomic write
              efactor[l][q][m][p] = eta;
            } else {
              if (phimin > 0.0) {
                double eta = eta_repulsive(phimin);
                //#pragma omp atomic write
                efactor[l][q][m][p] = eta;
              } else {
                double eta = eta_attractive(phimin);
                //#pragma omp atomic write
                efactor[l][q][m][p] = eta;
              }
            }
          }
        }
      }
    }

    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    BOOST_LOG_SEV(lg, info) << "Done computing enhancement factor...";
    BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();

    BOOST_LOG_SEV(lg, info)
        << "Computing neutral pairs : " << neutral_pairs.size();
    start = std::chrono::system_clock::now();

    for (unsigned int i = 0; i < neutral_pairs.size(); ++i) {
      short l = neutral_pairs[i].l_;
      short q = neutral_pairs[i].q_;
      short m = neutral_pairs[i].m_;
      short p = neutral_pairs[i].p_;
      efactor[l][q][m][p] = 1.0;
      efactor[m][p][l][q] = 1.0;
    }

    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    BOOST_LOG_SEV(lg, info) << "Done computing neutral pairs...";
    BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();
  }
  //////////////////////

  inline void compute_vdwenhancement_factor(double cutoff) {

    auto start = std::chrono::system_clock::now();

    BOOST_LOG_SEV(lg, info) << "Size for array efactor : " << efactor.size();

    BOOST_LOG_SEV(lg, info) << "Enhancement: expanding pair potentials";

// iterate in reduced_pairs
#pragma omp parallel for
    for (unsigned int irp = 0; irp < reduced_pairs.size(); ++irp) {
      // this is the index of reduced pairs, i.e., non repeated pair
      unsigned int index = reduced_pairs[irp].id_;

      // contact potential for this reduced pairhas the same ordering
      // of reduced pairs
      double contact_potential = contact_potentials[irp].potential_;

      // get repeated combinations (of particle_pairs)
      std::vector<long> reps = reduced_pairs[irp].repetitions_;

      // find if index is in barrier_potentials, i.e., if a potential
      // barrier exists
      ContactPotentialIterator barrier_it = std::find_if(
          barrier_potentials.begin(), barrier_potentials.end(), find_id(irp));

      // iterate in vector of repeated combinations
      for (unsigned int jrep = 0; jrep < reps.size(); ++jrep) {
        long rep_index = reps[jrep];
        short l = particle_pairs[rep_index].l_;
        short q = particle_pairs[rep_index].q_;
        short m = particle_pairs[rep_index].m_;
        short p = particle_pairs[rep_index].p_;

        double tr1 = rarray[l];  //(!notswapd ? rarray[l] : rarray[m]);
        // std::cerr << "\nn q is " << q;
        double tq1 = qarray[q];  //(!notswapd ? qarray[q] : qarray[p]);
        double potprefactor = tq1 * tq1 / tr1;
        double phimin = potprefactor * contact_potential;
        // std::cerr << "\n\n";

        cpotentials[m][p][l][q] = phimin;
        cpotentials[l][q][m][p] = phimin;
        // we have potential barrier
        // ** if barrier compute eta
        if (barrier_it != barrier_potentials.end()) {
          double phimax = potprefactor * (*barrier_it).potential_;
          // update potentials
          bpotentials[m][p][l][q] = phimax;
          bpotentials[l][q][m][p] = phimax;
          unsigned int idr = barrier_it - barrier_potentials.begin();
          rbarriers[m][p][l][q] = rbarrier_array[idr] * tr1;
          rbarriers[l][q][m][p] = rbarrier_array[idr] * tr1;
        }

        // std::cout << std::endl << l << '\t' << q << '\t' << m << '\t' << p <<
        // '\t' << efactor[l][q][m][p];
      }
    }

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    BOOST_LOG_SEV(lg, info) << "Done expanding pair potentials";
    BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();

    BOOST_LOG_SEV(lg, info) << "Computing Enhancement Factor";
    start = std::chrono::system_clock::now();

#pragma omp parallel for
    for (unsigned int l = 0; l < rarray_size; ++l) {
      // iterate in charges particle 1
      for (unsigned int q = 0; q < qarray_size; ++q) {
        // iterate in radii particle 2
        for (unsigned int m = 0; m < rarray_size; ++m) {
          // iterate in charges particle 2
          for (unsigned int p = 0; p < qarray_size; ++p) {
            double phimin = cpotentials[l][q][m][p]
                          + eint::potential_vdw(rarray[l]+rarray[m],
                                                rarray[l],
                                                rarray[m],
                                                cutoff);
            double phimax = bpotentials[l][q][m][p];

            // checks if a barrier exists over potential_threshold
            // FIXME change to 0
            // if (fabs(phimax) > potential_threshold) {
            if (fabs(phimax) != 0.0) {
              // double phimax = bpotentials[l][q][m][p];
              double eta = eta_barrier(phimin, phimax);

              // In the hybrid approach eta can be negative
              // in the case of phimin > phimax
              if (eta < 0.0) {
                BOOST_LOG_SEV(lg, warning) << "Negative enhancement factor";
                BOOST_LOG_SEV(lg, warning) << "phimin : " << phimin;
                BOOST_LOG_SEV(lg, warning) << "phimax : " << phimax;
                if (phimin > 0.0) {
                  eta = eta_repulsive(phimin);
                } else {
                  eta = eta_attractive(phimin);
                }
              }
              //#pragma omp atomic write
              efactor[l][q][m][p] = eta;
            }
            else {
              if (phimin > 0.0) {
                double eta = eta_repulsive(phimin);
                //#pragma omp atomic write
                efactor[l][q][m][p] = eta;
              } else {
                double eta = eta_attractive(phimin);
                //#pragma omp atomic write
                efactor[l][q][m][p] = eta;
              }
            }
          }
        }
      }
    }

    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    BOOST_LOG_SEV(lg, info) << "Done computing enhancement factor...";
    BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();

    BOOST_LOG_SEV(lg, info)
        << "Computing neutral pairs : " << neutral_pairs.size();
    start = std::chrono::system_clock::now();

    // for (unsigned int i = 0; i < neutral_pairs.size(); ++i) {
    //   short l = neutral_pairs[i].l_;
    //   short q = neutral_pairs[i].q_;
    //   short m = neutral_pairs[i].m_;
    //   short p = neutral_pairs[i].p_;
    //   efactor[l][q][m][p] = 1.0;
    //   efactor[m][p][l][q] = 1.0;
    // }

    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    BOOST_LOG_SEV(lg, info) << "Done computing neutral pairs...";
    BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();
  }
  //////////////////////

  inline void compute_enhancementfactor_frompairs() {
    auto start = std::chrono::system_clock::now();

    BOOST_LOG_SEV(lg, info) << "Size for array efactor : " << efactor.size();

    BOOST_LOG_SEV(lg, info) << "Enhancement: expanding pair potentials";

// iterate in reduced_pairs
#pragma omp parallel for
    for (unsigned int irp = 0; irp < reduced_pairs.size(); ++irp) {
      // this is the index of reduced pairs, i.e., non repeated pair
      unsigned int index = reduced_pairs[irp].id_;

      // contact potential for this reduced pairhas the same ordering
      // of reduced pairs
      double contact_potential = contact_potentials[irp].potential_;

      // get repeated combinations (of particle_pairs)
      std::vector<long> reps = reduced_pairs[irp].repetitions_;

      // find if index is in barrier_potentials, i.e., if a potential
      // barrier exists
      ContactPotentialIterator barrier_it = std::find_if(
          barrier_potentials.begin(), barrier_potentials.end(), find_id(irp));

      // iterate in vector of repeated combinations
      for (unsigned int jrep = 0; jrep < reps.size(); ++jrep) {
        long rep_index = reps[jrep];
        short l = particle_pairs[rep_index].l_;
        short q = particle_pairs[rep_index].q_;
        short m = particle_pairs[rep_index].m_;
        short p = particle_pairs[rep_index].p_;

        double tr1 = rarray[l];  //(!notswapd ? rarray[l] : rarray[m]);
        // std::cerr << "\nn q is " << q;
        double tq1 = qarray[q];  //(!notswapd ? qarray[q] : qarray[p]);
        double potprefactor = tq1 * tq1 / tr1;
        double phimin = potprefactor * contact_potential;
        // std::cerr << "\n\n";

        // we have potential barrier
        // ** if barrier compute eta
        if (barrier_it != barrier_potentials.end()) {
#pragma omp critical
          {
            // potential at contact
            PairElement bcpe1(m, p, l, q, phimin);
            bcpotentials_pe.push_back(bcpe1);
            PairElement bcpe2(l, q, m, p, phimin);
            bcpotentials_pe.push_back(bcpe2);

            double phimax = potprefactor * (*barrier_it).potential_;
            PairElement bpe1(m, p, l, q, phimax);
            bpotentials_pe.push_back(bpe1);
            PairElement bpe2(l, q, m, p, phimax);
            bpotentials_pe.push_back(bpe2);

            unsigned int idr = barrier_it - barrier_potentials.begin();
            double rbb = rbarrier_array[idr] * tr1;
            PairElement rpe1(m, p, l, q, rbb);
            rbarriers_pe.push_back(rpe1);
            PairElement rpe2(l, q, m, p, rbb);
            rbarriers_pe.push_back(rpe2);
          }
        } else {  // we dont have a barrier
#pragma omp critical
          {
            PairElement cpe1(m, p, l, q, phimin);
            cpotentials_pe.push_back(cpe1);
            PairElement cpe2(l, q, m, p, phimin);
            cpotentials_pe.push_back(cpe2);
          }
        }
        // std::cout << std::endl << l << '\t' << q << '\t' << m << '\t' << p <<
        // '\t' << efactor[l][q][m][p];
      }
    }

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    BOOST_LOG_SEV(lg, info) << "Done expanding pair potentials";
    BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();

    BOOST_LOG_SEV(lg, info) << "Computing Enhancement Factor";
    start = std::chrono::system_clock::now();

    // compute enhancement factor monotonic potential
#pragma omp parallel for  // ordered
    for (unsigned int imon = 0; imon < cpotentials_pe.size(); ++imon) {
      short l, q, m, p;
      double phimin;

      get_pe(cpotentials_pe[imon], l, q, m, p, phimin);

      if (phimin > 0.0) {
#pragma omp critical
        {
          double eta = eta_repulsive(phimin);
          PairElement eta_e(l, q, m, p, eta);
          efactor_pe.push_back(eta_e);
        }
      } else {
#pragma omp critical
        {
          double eta = eta_attractive(phimin);
          PairElement eta_e(l, q, m, p, eta);
          efactor_pe.push_back(eta_e);
        }
      }
    }

    // compute enhancement factor barrier potential
#pragma omp parallel for  // ordered
    for (unsigned int ibar = 0; ibar < rbarriers_pe.size(); ++ibar) {
      short l, q, m, p;
      short lx, qx, mx, px;
      double phimin, phimax;

      get_pe(bcpotentials_pe[ibar], l, q, m, p, phimin);
      get_pe(bpotentials_pe[ibar], lx, qx, mx, px, phimax);

      // checks if a barrier exists over potential_threshold
      // FIXME change to 0
      // if (fabs(phimax) > potential_threshold) {
      if (fabs(phimax) != 0.0) {
        // double phimax = bpotentials[l][q][m][p];
        double eta = eta_barrier(phimin, phimax);

        // In the hybrid approach eta can be negative
        // in the case of phimin > phimax
        if (eta < 0.0) {
          BOOST_LOG_SEV(lg, warning) << "Negative enhancement factor";
          BOOST_LOG_SEV(lg, warning) << "phimin : " << phimin;
          BOOST_LOG_SEV(lg, warning) << "phimax : " << phimax;
          if (phimin > 0.0) {
            eta = eta_repulsive(phimin);
          } else {
            eta = eta_attractive(phimin);
          }
        }
        //#pragma omp atomic write

#pragma omp critical
        {
          PairElement eta_e(l, q, m, p, eta);
          efactor_pe.push_back(eta_e);
        }
      }
    }

    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    BOOST_LOG_SEV(lg, info) << "Done computing enhancement factor...";
    BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();

    BOOST_LOG_SEV(lg, info)
        << "Computing neutral pairs : " << neutral_pairs.size();
    start = std::chrono::system_clock::now();

    for (unsigned int i = 0; i < neutral_pairs.size(); ++i) {
      short l = neutral_pairs[i].l_;
      short q = neutral_pairs[i].q_;
      short m = neutral_pairs[i].m_;
      short p = neutral_pairs[i].p_;
      PairElement eta_1(l, q, m, p, 1.0);
      efactor_pe.push_back(eta_1);
      PairElement eta_2(m, p, l, q, 1.0);
      efactor_pe.push_back(eta_2);
    }

    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    BOOST_LOG_SEV(lg, info) << "Done computing neutral pairs...";
    BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();

    BOOST_LOG_SEV(lg, info)
        << "Enhancement factor to arrays : " << efactor_pe.size();
    start = std::chrono::system_clock::now();

    bsgrid2d = {{static_cast<long int>(efactor_pe.size()), 4}};
    efindices.resize(bsgrid2d);
    daefactor.resize(efactor_pe.size());

    //#pragma omp parallel for ordered
    for (unsigned int ieta = 0; ieta < efactor_pe.size(); ++ieta) {
      short l, q, m, p;
      double eta;
      get_pe(efactor_pe[ieta], l, q, m, p, eta);
      efindices[ieta][0] = l;
      efindices[ieta][1] = q;
      efindices[ieta][2] = m;
      efindices[ieta][3] = p;
      daefactor[ieta] = eta;
      // std::cerr << "\n" << l << '\t'  << q << '\t'  << m << '\t'  << p <<
      // '\t'  << eta;
    }

    // Contact monotonic potentials
    bsgrid2d = {{static_cast<long int>(cpotentials_pe.size()), 4}};
    cpindices.resize(bsgrid2d);
    dacpotentials.resize(cpotentials_pe.size());

    for (unsigned int icp = 0; icp < cpotentials_pe.size(); ++icp) {
      short l, q, m, p;
      double cpot;
      get_pe(cpotentials_pe[icp], l, q, m, p, cpot);
      cpindices[icp][0] = l;
      cpindices[icp][1] = q;
      cpindices[icp][2] = m;
      cpindices[icp][3] = p;
      dacpotentials[icp] = cpot;
      // std::cerr << "\n" << l << '\t'  << q << '\t'  << m << '\t'  << p <<
      // '\t'  << eta;
    }

    // Barrier potentials
    bsgrid2d = {{static_cast<long int>(rbarriers_pe.size()), 4}};
    bpindices.resize(bsgrid2d);
    darbarriers.resize(rbarriers_pe.size());
    dabpotentials.resize(rbarriers_pe.size());
    dabcpotentials.resize(rbarriers_pe.size());

    for (unsigned int ibar = 0; ibar < rbarriers_pe.size(); ++ibar) {
      short l, q, m, p;
      double rbar;
      get_pe(rbarriers_pe[ibar], l, q, m, p, rbar);
      bpindices[ibar][0] = l;
      bpindices[ibar][1] = q;
      bpindices[ibar][2] = m;
      bpindices[ibar][3] = p;
      darbarriers[ibar] = rbar;
      dabpotentials[ibar] = bpotentials_pe[ibar].value_;
      dabcpotentials[ibar] = bcpotentials_pe[ibar].value_;
      // std::cerr << "\n" << l << '\t'  << q << '\t'  << m << '\t'  << p <<
      // '\t'  << eta;
    }

    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    BOOST_LOG_SEV(lg, info) << "Done enhancement factor to arrays...";
    BOOST_LOG_SEV(lg, info) << "Elapsed time : " << elapsed_seconds.count();
  }

  inline void get_efindices(darray &daefactor_,
                            boost_short_array2d &efindices_) {
    // get number of elements
    long int efsize = static_cast<long int>(efindices.shape()[0]);

    // resize daefactor and copy elements
    daefactor_.resize(efsize);
    std::copy(std::begin(daefactor), std::end(daefactor),
              std::begin(daefactor_));

    // resize indices and copy elements
    bsgrid2d = {{efsize, 4}};
    efindices_.resize(bsgrid2d);
    std::copy(efindices.begin(), efindices.end(), efindices_.begin());
  }

  inline void get_cpindices(darray &dacpotentials_,
                            boost_short_array2d &cpindices_) {
    // get number of elements
    long int cpsize = static_cast<long int>(cpindices.shape()[0]);

    // resize vectors and copy elements
    dacpotentials_.resize(cpsize);
    std::copy(std::begin(dacpotentials), std::end(dacpotentials),
              std::begin(dacpotentials_));

    // resize indices and copy elements
    bsgrid2d = {{cpsize, 4}};
    cpindices_.resize(bsgrid2d);
    std::copy(cpindices.begin(), cpindices.end(), cpindices_.begin());
  }

  inline void get_bpindices(darray &dabcpotentials_, darray &dabpotentials_,
                            darray &darbarriers_,
                            boost_short_array2d &bpindices_) {
    // get number of elements
    long int bpsize = static_cast<long int>(bpindices.shape()[0]);

    // resize vectors and copy elements
    dabcpotentials_.resize(bpsize);
    std::copy(std::begin(dabcpotentials), std::end(dabcpotentials),
              std::begin(dabcpotentials_));

    dabpotentials_.resize(bpsize);
    std::copy(std::begin(dabpotentials), std::end(dabpotentials),
              std::begin(dabpotentials_));

    darbarriers_.resize(bpsize);
    std::copy(std::begin(darbarriers), std::end(darbarriers),
              std::begin(darbarriers_));

    // resize indices and copy elements
    bsgrid2d = {{bpsize, 4}};
    bpindices_.resize(bsgrid2d);
    std::copy(bpindices.begin(), bpindices.end(), bpindices_.begin());
  }

  //////////////////////
  darray rarray;
  darray qarray;

  // no ref
  boost_array4d_ref efactor;
  boost_array4d_ref cpotentials;
  boost_array4d_ref bpotentials;
  boost_array4d_ref rbarriers;

  boost_array4d dummy;

  boost_short_array4d nterms4d;

  darray daefactor;
  darray dacpotentials;
  darray dabpotentials;
  darray dabcpotentials;
  darray darbarriers;

  double eps;
  double AH;

  src::severity_logger<severity_level> lg;

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
  std::vector<ContactPotential> error_potentials;
  std::vector<ContactPotential> barrier_potentials;

  std::vector<PairElement> efactor_pe;
  std::vector<PairElement> cpotentials_pe;  //!< contact potential no barrier
  std::vector<PairElement>
      bcpotentials_pe;  //!< contact potential when barrier exists
  std::vector<PairElement> bpotentials_pe;  //!< barrier potential
  std::vector<PairElement> rbarriers_pe;    //!< barrier location

  boost_short_array2d efindices;  //!< Array for enhancement factor indices.
  boost_short_array2d cpindices;  //!< Array for contact potential indices.
  boost_short_array2d
      bpindices;  //!< Array for rbarrier and barrier potential indices.

  bshortgrid2d bsgrid2d;  //!< Grid for indices

  double potential_threshold;
  double temperature = 300.0;
  std::vector<double> rbarrier_array;
};
}  // namespace enhancement

#endif  // ENHANCEMENT_H
