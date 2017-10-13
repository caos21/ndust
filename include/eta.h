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
  unsigned int l_;            //!< Index l of current volume section
  unsigned int q_;            //!< Index q of current charge section
  unsigned int m_;            //!< Index m of coalescing volume
  unsigned int p_;            //!< Index p of coalescing charge
  unsigned int n_;            //!< Index n of coalescing volume
  unsigned int r_;            //!< Index r of coalescing volume
  double eta_;                //!< Value for \f$\eta_{mp,nr}(v_l, Q_q)\f$

  friend void fill_eta(EtaCreationFactor &ecf,
                       const unsigned int id,
                       const unsigned int l,
                       const unsigned int q,
                       const unsigned int m,
                       const unsigned int p,
                       const unsigned int n,
                       const unsigned int r,
                       const double eta);

  friend void etafactor_tostream(const EtaCreationFactor &ecf,
                                 std::fstream& etafactor_file);

  friend void etafactor_fromstream(EtaCreationFactor &ecf,
                                   std::istringstream &ss);
};

inline void fill_eta(EtaCreationFactor &ecf,
                     const unsigned int id,
                     const unsigned int l,
                     const unsigned int q,
                     const unsigned int m,
                     const unsigned int p,
                     const unsigned int n,
                     const unsigned int r,
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
  unsigned int l_;            //!< Index l of current volume section
  unsigned int q_;            //!< Index q of current charge section
  unsigned int m_;            //!< Index m of coalescing volume
  unsigned int p_;            //!< Index p of coalescing charge
  unsigned int n_;            //!< Index n of coalescing volume
  unsigned int r_;            //!< Index r of coalescing volume
  double eta_;                //!< Value for \f$\eta_{mp,nr}(v_l, Q_q)\f$

  ss >> sid_
     >> sl_
     >> sq_
     >> sm_
     >> sp_
     >> sn_
     >> sr_
     >> seta_;

  id_  = std::stoul(sid_ );
  l_   = std::stoul(sl_  );
  q_   = std::stoul(sq_  );
  m_   = std::stoul(sm_  );
  p_   = std::stoul(sp_  );
  n_   = std::stoul(sn_  );
  r_   = std::stoul(sr_  );
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
  unsigned int l_;            //!< Index l of current volume section
  unsigned int q_;            //!< Index q of current charge section
  unsigned int m_;            //!< Index m of coalescing volume
  unsigned int p_;            //!< Index p of coalescing charge
  double death_;              //!< Value for death term

  friend void fill_death(DeathFactor &dth,
                       unsigned int id,
                       unsigned int l,
                       unsigned int q,
                       unsigned int m,
                       unsigned int p,
                       double death);

  friend void deathfactor_tostream(const DeathFactor &dth,
                                   std::fstream& deathfactor_file);

  friend void deathfactor_fromstream(DeathFactor &df,
                                     std::istringstream &ss);
};

inline
void fill_death(DeathFactor &dth,
                unsigned int id,
                unsigned int l,
                unsigned int q,
                unsigned int m,
                unsigned int p,
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
  unsigned int l_;            //!< Index l of current volume section
  unsigned int q_;            //!< Index q of current charge section
  unsigned int m_;            //!< Index m of coalescing volume
  unsigned int p_;            //!< Index p of coalescing charge
  double death_;              //!< Value for \f$\eta_{mp,nr}(v_l, Q_q)\f$

  ss >> sid_
     >> sl_
     >> sq_
     >> sm_
     >> sp_
     >> sdeath_;

  id_  = std::stoul(sid_ );
  l_   = std::stoul(sl_  );
  q_   = std::stoul(sq_  );
  m_   = std::stoul(sm_  );
  p_   = std::stoul(sp_  );
  death_ = std::stod(sdeath_);

  fill_death(df, id_, l_, q_, m_, p_, death_);
}
#endif
