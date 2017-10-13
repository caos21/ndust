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

#ifndef IO_H
#define IO_H

#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <exception>
#include <string>
#include <boost/variant.hpp>
#include <cstdlib>

// ToDo:
//   * Documentation
//   * Fix for string parameters
//   * Fix for integers

// Get environment variable
inline
std::string get_environment_variable(std::string const& envar) {
  char *value = getenv(envar.c_str());
  return (value == NULL? std::string(""): std::string(value));
}


// Boost variant type, handles double unsigned bool and string
typedef boost::variant<double, unsigned int, bool, std::string> boost_variant;

// void print_header(std::ostream& stream);
//
// void print_parameters(std::ostream& stream,
//                       const std::map<std::string, boost_variant > mappar);
//
// char* get_cmd_option(char** begin, char** end, const std::string & option);
//
// bool cmd_option_exists(char** begin, char** end, const std::string& option);

// template <class T>
// void param_input(const char type, T *par_Value, const T def_Value,
//                  const std::string & option_Value, int argc, char** argv);
// void show_help(int argc, char** argv);

inline
void print_header(std::ostream& stream) {
  stream << std::endl;
  stream << std::endl << "####################################################";
  stream << std::endl << "######                 CRat                  #######";
  stream << std::endl << "####################################################";
  stream << std::endl;
  stream << std::endl;
}

inline
void print_parameters(std::ostream& stream, const std::map<std::string,
                      boost_variant > mappar) {
  std::map<std::string, boost_variant >::const_iterator it;

  for (it = mappar.begin(); it != mappar.end(); ++it) {
      stream << std::endl << "#<] " << it->first << " = " << it->second;
//              << std::endl;
  }
  stream << std::endl;
}

inline
char* get_cmd_option(char** begin, char** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

inline
bool cmd_option_exists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}

template <class T>
void param_input(const char type, T *par_Value, const T def_Value,
                 const std::string & option_Value, int argc, char** argv)
{
  // Store input option
  char *input_opt;

  // Check if option  is called
  input_opt = get_cmd_option(argv, argv + argc, option_Value);

  // Add exceptions bad_cast

  if(input_opt)
  {
    switch(type)
    {
//       // Integer value for parameter
//       case 'i': case 'I':
//         *par_Value = std::atoi(input_opt);
//         break;

      // Unsigned long
      case 'u': case 'U':
        *par_Value = std::strtoul(input_opt, 0, 10);
        break;

      // Double value
      case 'd': case 'D': case 'f': case 'F':
        *par_Value = std::atof(input_opt);
        break;
/*
      case 's': case 'S':
//         std::string aux(input_opt);
// //         *par_Value=*input_opt;
//         par_Value(input_opt);
//         std::string a;
//         a.
        par_Value->assign(input_opt, strlen(input_opt));
        break;*/

      default:
        std::cout << std::endl << "!! Error: bad type input for param_input. "
                  << std::endl;
        exit(-1);

    };
    std::cout << std::endl << "!! Using " << option_Value << " = "
              << *par_Value << std::endl;
  }
  else
  {
    *par_Value = def_Value;
    std::cout << std::endl << "!! Using default value "
              << option_Value << " = " << *par_Value << std::endl;
  }
}

inline
void show_help(int argc, char** argv)
{
  if(cmd_option_exists(argv, argv+argc, "-h"))
  {
    std::cout << std::endl << "Options:" << std::endl;
    std::cout << '\t' << "-d Directory for input/output files" << std::endl;
    std::cout << '\t' << "-o Prefix for input/output files" << std::endl;
    std::cout << '\t' << "-t Number of OpenMP threads" << std::endl;
    std::cout << std::endl;
    exit(0);
  }
}


#endif // IO_H
