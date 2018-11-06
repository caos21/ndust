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

#ifndef NDVERSION_H
#define NDVERSION_H

#include <iostream>

#include "log.h"

#include "ndust_version.h"

namespace ndversion {

  /**
  * class NDVersion
  *
  * Represents NDust version from git
  *
  */
  class NDVersion {
  public:
    //
    //! Constructor for NDVersion
    /*!
      @param  lg_ Logger instance.
    */
  NDVersion(src::severity_logger< severity_level > lg_): lg(lg_) {
    }

    void print() {
      std::cout << "Version: "<< g_NDUST_VERSION << std::endl;
      std::cout << "Hash: " << g_NDUST_HASH << std::endl;
      std::cout << "State: "<< g_NDUST_STATE << std::endl;
      std::cout << "Origin: "<< g_NDUST_BRANCH << std::endl;
    }

    void log() {
      BOOST_LOG_SEV(lg, info) << "Version: "<< g_NDUST_VERSION;
      BOOST_LOG_SEV(lg, info) << "Hash: " << g_NDUST_HASH;
      BOOST_LOG_SEV(lg, info) << "State: " << g_NDUST_STATE;
      BOOST_LOG_SEV(lg, info) << "Origin: "<< g_NDUST_BRANCH;
    }
  private:
    src::severity_logger< severity_level > lg; //!< Logger instance.
  };
}

#endif// NDVERSION_H
