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

#ifndef DESCRIPTION_H
#define DESCRIPTION_H

/**
  * struct Description
  *
  * Stores the data description
  *
  */
struct description {
  std::string description;               //!< Run description
  std::string timestamp;                 //!< Timestamp of creation
  std::string sysname;                   //!< OS name
  std::string nodename;                  //!< Hostname
  std::string release;                   //!< kernel version
  std::string version;                   //!< OS version
  std::string machine;                   //!< Architecture
};

#endif
