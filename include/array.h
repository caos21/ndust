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

#ifndef ARRAY_H
#define ARRAY_H

#include <valarray>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cmath>

#ifdef USING_EIGEN

#include <Eigen/Eigen>
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixDD;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorDD;

#endif

// Boost array
#include <boost/multi_array.hpp>
#include <boost/math/special_functions/cbrt.hpp>

typedef boost::multi_array<double, 1> boost_array1d;

typedef boost::multi_array<double, 2> boost_array2d;
typedef boost::multi_array<double, 3> boost_array3d;
typedef boost::multi_array<double, 4> boost_array4d;

typedef boost::multi_array_ref<double, 4> boost_array4d_ref;
typedef boost::const_multi_array_ref<double, 4> const_boost_array4d_ref;

typedef boost::array<boost_array2d::index, 2> bgrid2d;
typedef boost::array<boost_array3d::index, 3> bgrid3d;
typedef boost::array<boost_array4d::index, 4> bgrid4d;

typedef boost::multi_array<short, 2> boost_short_array2d;
typedef boost::array<boost_short_array2d::index, 2> bshortgrid2d;
typedef boost::multi_array_ref<short, 2> boost_shortarray2d_ref;

// u integer type
typedef boost::multi_array<unsigned int, 2> boost_uint_array2d;
typedef boost::array<boost_uint_array2d::index, 2> buigrid2d;

typedef boost::multi_array<short unsigned int, 4> boost_short_array4d;

// Define template array = valarray
template<typename Type>
using array = std::valarray< Type >;

// Define darray as an array of doubles
typedef array<double> darray;
// typedef std::valarray<double> darray;
// typedef std::vector<double> darray;

// WARNING DANGER FIXME boost::multi_array<TYPE, 2>
// boost::multi_array for an array of two indices
template<typename Type>
using array2d = boost::multi_array<double, 2>;

// defines the 2d grid
template<typename Type>
using grid2d = boost::array<array2d<double>::index, 2>;

// Define darray2d as an array2d of doubles
typedef array2d<double> darray2d;

// defines grid2d for doubles
typedef grid2d<double> dgrid2d;

// inspired by http://stackoverflow.com/questions/27028226/python-linspace-in-c
template<typename Type>
darray linear(Type start_in, Type end_in, unsigned int num_in) {
  // Cast
  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  // calculate delta
  double delta = (end - start) / (num - 1.0);

  // create array
  darray linear(num);

  // populate array
  for(unsigned int i=0; i < num; ++i) {
    linear[i] = start + delta * i;
  }
  linear[num-1] = end;

  return linear;
}

inline
int find_index(const darray array, double value) {
  int index = std::distance(std::begin(array),
			    std::find(std::begin(array),
				      std::end(array), value));
  if (index == array.size()) {
    return -1;
  }
  else {
    return index;
  }
}

inline
double mean(const darray array) {
  return array.sum()/static_cast<double>(array.size());
}

// return the maximum in absolute value
template <class T> const T& maxabs(const T& a, const T& b) {
  double abs_a = abs(a);
  double abs_b = abs(b);
  return (abs_a<abs_b)?abs_b:abs_a;
}

// return the minimum in absolute value
template <class T> const T& minabs(const T& a, const T& b) {
  double abs_a = abs(a);
  double abs_b = abs(b);
  return !(abs_b<abs_a)?abs_a:abs_b;
}


template <class ForwardIterator>
ForwardIterator max_abselem(ForwardIterator first, ForwardIterator last) {
  if (first==last) return last;
  ForwardIterator largest = first;

  while (++first!=last)
    if (abs(*largest)<abs(*first))
      largest=first;

  return largest;
}

template <class ForwardIterator>
ForwardIterator min_abselem(ForwardIterator first, ForwardIterator last) {
  if (first==last) return last;
  ForwardIterator smallest = first;

  while (++first!=last)
    if (!(abs(*smallest)<abs(*first)))
      smallest=first;
  return smallest;
}

// Computes the average of an array 1/2(a_i + a_i+1)
inline darray average(darray array) {
  // create valarray

  unsigned int size = array.size();

  darray avg(size-1);

// WARNING possible implementation using slices start, numelem, stride
//   darray hi = array[std::slice(1, size-1, 1)];
//   darray low = array[std::slice(0, size-1, 1)];
//   avg = ( hi + low );

  for(unsigned int i=0; i < size-1; ++i) {
    avg[i] = array[i] + array[i+1];
  }

  return 0.5*avg;
}
// def average(nparray):
//     """ Computes the average of an array 1/2(a_i + a_i+1)
//     """
//     return 0.5 * (nparray[1:] + nparray[:-1])

template < class ArrType >void PrintArray(ArrType varray, unsigned int size) {
  //
  std::cout << std::endl;
  for (unsigned int i = 0; i < size; i++) {
    std::cout << varray[i] << std::endl;
  }
  std::cout << "# Size" << size << std::endl;
}

template < class ArrType >void PrintArray(ArrType varray) {
  //
  unsigned int size = varray.size();
  std::cout << std::endl;
  for (unsigned int i = 0; i < size; i++) {
    std::cout << varray[i] << std::endl;
  }
}

template < class ArrType >void WriteArray(const ArrType varray,
                                           unsigned int size,
                                           std::fstream* file_stream,
                                           const std::string message = "") {
  //
  if (!message.empty())
    *file_stream << message << std::endl;
  for (unsigned int i = 0; i < size; i++) {
    *file_stream << varray[i] << std::endl;
  }
}

template < class ArrType, class Type >
ArrType ReadArray(std::fstream* file_stream) {
  //
//   ArrType array;
  std::vector<Type> vecarray;
  Type val;
  unsigned int n = 0;
  while (*file_stream >> val) {
    vecarray.push_back(val);
    ++n;
  }
  ArrType array(n);
  for(unsigned int i=0; i<n; ++i) {
    array[i] = vecarray[i];
  }

  return array;
}

template < class Type >
void PrintArray2d(array2d<Type> arr2d) {
  std::cout << std::endl << "[ii] === PrintArray2d ===" << std::endl;
  for(size_t x = 0; x < arr2d.shape()[0]; ++x) {
      std::cout << "[ii] ";
      for(size_t y = 0; y < arr2d.shape()[1]; ++y) {
        std::cout << arr2d[x][y] << '\t';
      }
      std::cout << std::endl;
  }
  std::cout << std::endl;
}

inline void CopyArray2d(const darray2d& source, darray2d& dest) {

  dest.resize(boost::extents[source.shape()[0]][source.shape()[1]]);
  dest = source;
}

template < class Type >
inline void CopyArray2d(const Type& source, Type& dest) {
  dest.resize(boost::extents[source.shape()[0]][source.shape()[1]]);
  dest = source;
}

template < class Type >
inline void CopyArray2d(const array2d<double>& source, array2d<double>& dest) {
  // get shape source array and reshape dest
//   std::vector<size_t> grid;
//   const size_t* shape = source.shape();
//   grid.assign(shape, shape+source.num_dimensions());
//   dest.resize(grid);
//   for(unsigned int i=0; i<source.shape()[0]; ++i) {
//     for(unsigned int j=0; j<source.shape()[1]; ++j) {
//       dest[i][j] = source[i][j];
//     }
//   }
  dest.resize(boost::extents[source.shape()[0]][source.shape()[1]]);
  dest = source;
}

// WARNING DELETE
// inline darray2d CopyArray2d(const darray2d source) {
//   darray2d dest = new darray2d(source);
// //   return &dest;
// }

template < class Type >
array2d<Type> ReadArray2d(std::fstream* file_stream, std::string& comments) {
  //
  // WARNING FIXME use one vector or array2d
  std::vector<Type> xvecarray;
  std::vector<Type> yvecarray;
  Type valx, valy;
  unsigned int n = 0;
  unsigned int nlines = 0;

  std::stringstream sscomm;
  // inspiration http://stackoverflow.com/questions/576677/how-do-i-skip-reading-a-line-in-a-file-in-c
  std::string sline;
  while (std::getline(*file_stream, sline)) {
    // Trims white space at the beginning of the string
    sline.erase(sline.begin(),
                std::find_if(sline.begin(), sline.end(),
                             std::not1(std::ptr_fun<int, int>(isspace))));

    if(sline[0] == '#') {
      sscomm << "[+" << nlines << "] " << sline << std::endl;
    }
    else {
      std::stringstream(sline) >> valx >> valy;
      xvecarray.push_back(valx);
      yvecarray.push_back(valy);
      ++n;
    }
    ++nlines;
  }

  comments = sscomm.str();
  std::cout  << std::endl << "[ii] " << nlines << " lines readed";

  // WARNING FIXME template see comments above
  grid2d<Type> table_t = {{n, 2}};
  array2d<Type> table(table_t);

  for(unsigned int i=0; i<n; ++i) {
    table[i][0] = xvecarray[i];
    table[i][1] = yvecarray[i];
  }

  return table;
}

template < class ArrType1, class ArrType2 >
void Write2Arrays(const ArrType1 array1, const ArrType2 array2,
                  unsigned int size, std::fstream* file_stream,
                  const std::string message = "") {
  //
  if (!message.empty())
    *file_stream << message << std::endl;
  for (unsigned int i = 0; i < size; i++) {
    *file_stream << array1[i] << '\t' << array2[i] << std::endl;
  }
}

// WARNING FIXME template for handle files and streams ???
//   Write2Arrays<darray, darray>(x, y, 200, &std::cout, "\n#test");
template < class ArrType1, class ArrType2 >
void Write2Arrays(const ArrType1 array1, const ArrType2 array2,
                  unsigned int size, std::ostream* out_stream,
                  const std::string message = "") {
  //
  if (!message.empty())
    *out_stream << message << std::endl;
  for (unsigned int i = 0; i < size; i++) {
    *out_stream << array1[i] << '\t' << array2[i] << std::endl;
  }
}

template < class ArrType1, class ArrType2 > double Trapeze(ArrType1 xarray,
                                          ArrType2 yarray, unsigned int size) {
  double x0 = xarray[0];
  double xn = xarray[size-1];

  double sum = 0.0;

  for (unsigned int i = 0; i < size-1; i++) {
    sum += (xarray[i+1]-xarray[i]) * (yarray[i+1]+yarray[i]);
  }
  return sum;
}

inline void mod(std::valarray<double> &da) {
  da *= 2.0;
}

inline
bool computeProfiles(std::vector<darray> newprof,
                     std::vector<darray> oldprof,
                     unsigned int npoints,
                     darray& meanvec,
                     double profile_rerr) {
  meanvec = 0.0;
  std::vector<darray>::iterator np;
  std::vector<darray>::iterator op;
  for(unsigned int i=0; i<npoints; i++) {
    unsigned int nd = 0;
    for(np=newprof.begin(), op=oldprof.begin();
        np!=newprof.end(); ++np, ++op) {
      double npval = (*np)[i];
      double opval = (*op)[i];
      meanvec[nd] += npval;
      if(abs(npval - opval) > std::min(npval, opval) * profile_rerr) {
//         std::cout << std::endl << "Continue "
//                   << abs(p1val - p2val) / std::min(p1val, p2val);
        std::cout << std::endl << "[ii] Profiles differ " << std::endl;
        return false;
      }
      ++nd;
    }
  }
  meanvec /= static_cast<double>(npoints);
  std::cout << std::endl << "[ii] Profile comparison successful " << std::endl;
  return true;
}

inline
void normalize_profile(const darray profile, darray& normalized) {
  normalized = profile / profile.max();
  std::valarray<bool> nans = (normalized != normalized);
  normalized[nans] = 0.0;
}

template < typename Type >
inline Type vabs(Type elem) {
  return abs(elem);
}


inline
bool almost_equal(double a, double b, double tol=1e-5) {
//   if(abs(a-b) < std::min(abs(a), abs(b)) * 0.01) return true;
  if(abs(a-b) < tol) return true;
  else return false;
}


inline
bool is_close(double a, double b, double rel_tol=1e-9, double abs_tol=0.0) {
  if (std::abs(a-b) <= std::max(rel_tol * std::max(std::abs(a), std::abs(b)), abs_tol)) {
    return true;
  }
  else {
    return false;
  }	    
}

#endif// ARRAY_H
