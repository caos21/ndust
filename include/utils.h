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

#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <utility>
#include <string>
#include <functional>
#include <algorithm>
#include <parallel/algorithm>
#include <cstdlib>

#include <omp.h>

// multiprecision
#include <bits/stdc++.h>
#include <boost/multiprecision/cpp_int.hpp>
using boost::multiprecision::cpp_int;
#include <boost/multiprecision/cpp_dec_float.hpp>
using boost::multiprecision::cpp_dec_float_50;
using boost::multiprecision::cpp_dec_float_100;

// from https://www.geeksforgeeks.org/factorial-large-number-using-boost-multiprecision-library/
inline
cpp_int Factorial(int number) {
  cpp_int num = 1;
  for (int i = 1; i <= number; i++)
    num = num * i;
  return num;
}

#include <omp.h>


namespace utils {

  template <typename T>
  T kron_delta(const T i, const T j) {
    return (i == j ? T(1) : T(0));
  }

  template <typename T>
    T absolute_error(T true_value, T expected_value){
    return expected_value - true_value;
  }
  
  template <typename T>
    T relative_error(T true_value, T expected_value){
    return fabs((expected_value - true_value)/true_value);
  }

  template <typename T>
    T max_relative_error(T a, T b){
    return std::max(relative_error(a, b), relative_error(b, a));
  }

  template <typename T>
    T max_pct_error(T a, T b){
    return 100.0 * max_relative_error(a, b);
  }
 
}

#endif//UTILS_H




