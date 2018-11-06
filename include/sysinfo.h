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

#ifndef SYSINFO_H
#define SYSINFO_H

#include <iostream>
#include <sys/utsname.h>

#include "log.h"

namespace sysinfo {

  /**
  * class SysInfo
  *
  * Represents system information strings
  *
  */
  class SysInfo {
  public:
    //
    //! Constructor for SysInfo
    /*!
      @param  lg_ Logger instance.
    */
  SysInfo(src::severity_logger< severity_level > lg_): lg(lg_) {
      uname(&suname);
    }

    void print(){
      std::cout << "System Name: "<< suname.sysname << std::endl;
      std::cout << "Hostname: "<< suname.nodename << std::endl;
      std::cout << "Kernel version: "<< suname.release << std::endl;
      std::cout << "Kernel build timestamp: "<< suname.version << std::endl;
      std::cout << "Machine architecture: "<< suname.machine << std::endl;
      std::cout << "Domain name: "<< suname.domainname << std::endl;
    }

    void log(){
      BOOST_LOG_SEV(lg, info) << "System Name: "<< suname.sysname;
      BOOST_LOG_SEV(lg, info) << "Hostname: "<< suname.nodename;
      BOOST_LOG_SEV(lg, info) << "Kernel version: "<< suname.release;
      BOOST_LOG_SEV(lg, info) << "Kernel build timestamp: "<< suname.version;
      BOOST_LOG_SEV(lg, info) << "Machine architecture: "<< suname.machine;
      BOOST_LOG_SEV(lg, info) << "Domain name: "<< suname.domainname;
    }
  private:
    src::severity_logger< severity_level > lg; //!< Logger instance.
    struct utsname suname;    
  };
}

#endif// SYSINFO_H
