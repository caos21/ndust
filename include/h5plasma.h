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

#ifndef H5PLASMA_H
#define H5PLASMA_H

#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <exception>
#include <string>
#include <cstdlib>
#include <cassert>
#include <cstring>

#include <H5Cpp.h>

#include "array.h"

#define H5_DSET_ERR -5

// #define VERBOSE

// TODO add namespace
// TODO overloading and template specialization


inline
int create_group_hdf5(H5::H5File file,
                      std::string sgroup) {
  
  try {
    // Turn off the auto-printing
    H5::Exception::dontPrint();

    H5::Group group;
    // attempt to open group
    try {
      H5::Exception::dontPrint();
      group = file.openGroup(sgroup);
    }
    catch( H5::Exception gerror ) {
//       gerror.printError();
      // if group does not exist, create it
      try {
        H5::Exception::dontPrint();
        group = file.createGroup(sgroup);
      }
      catch( H5::Exception gerror2 ) {
        gerror2.printError();
        std::cout << std::endl << "[ee] h5 I/O error in group "
                  << sgroup << ". Terminate.\n";
        std::terminate();
      }
    }
  }
  // catch failure caused by the H5File operations
  catch( H5::FileIException error )
  {
    error.printError();
    return -1;
  }
  return 0;
}

inline
int write_dataset2d_hdf5(H5::H5File file,
                         std::string sgroup,
                         std::string sdataset,
                         boost_array2d array2D,
                         unsigned int nrow,
                         unsigned int ncol) {

  double* varray = new double[nrow*ncol];
  for(unsigned int i = 0; i<nrow; ++i) {
//     std::cout << std::endl;
    for(unsigned int j = 0; j<ncol; ++j) {
        varray[i+nrow*j] = array2D[i][j];
//           std::cout << varray[i+nrow*j] << '\t';
    }
  }
//   std::cout << std::endl << nrow << '\t' << ncol << std::endl;
  try
  {
    // Turn off the auto-printing
    H5::Exception::dontPrint();

    H5::Group group;
    // attempt to open group
    try {
      H5::Exception::dontPrint();
      group = file.openGroup(sgroup);
    }
    catch( H5::Exception gerror ) {
//       gerror.printError();
      // if group does not exist, create it
      try {
        H5::Exception::dontPrint();
        group = file.createGroup(sgroup);
      }
      catch( H5::Exception gerror2 ) {
        gerror2.printError();
        std::cout << std::endl << "[ee] h5 I/O error in group "
                  << sgroup << ". Terminate.\n";
        std::terminate();
      }
    }

    H5::DataSet dataset;
    // attempt to open group
    try {
      H5::Exception::dontPrint();
      dataset = group.openDataSet(sdataset);
    }
    catch( H5::Exception aerror )
    {
//       aerror.printError();
      try {
        H5::Exception::dontPrint();
        H5::DataType type(H5::PredType::NATIVE_DOUBLE);
        // dataset dimensions
        hsize_t dimsf[2];
        dimsf[0] = ncol;
        dimsf[1] = nrow;
        H5::DataSpace idset(2, dimsf);
        dataset = group.createDataSet(sdataset, type, idset);
        type.close();
        idset.close();
      }
      catch( H5::Exception aerror2 )
      {
        aerror2.printError();
        std::cout << std::endl << "[ee] h5 I/O error in dataset "
                  << sdataset << ". Terminate.\n";
        std::terminate();
      }
    }
//     // dataset dimensions
//     hsize_t dimsf[2];
//     dimsf[0] = ncol;
//     dimsf[1] = nrow;
//     H5::DataSpace dataspace(2, dimsf);

    H5::DataType datatype(H5::PredType::NATIVE_DOUBLE);

    dataset.write(varray, datatype);

    delete [] varray;

    // WARNING FIXME dataspace and group not closed
    datatype.close();
    dataset.close();
//     dataspace.close();
   }  // end of try block
   // catch failure caused by the H5File operations
   catch( H5::FileIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataSet operations
   catch( H5::DataSetIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataSpace operations
   catch( H5::DataSpaceIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataSpace operations
   catch( H5::DataTypeIException error )
   {
      error.printError();
      return -1;
   }
   return 0;
}


inline
int read_dataset2d_hdf5(H5::H5File file,
                        std::string sgroup,
                        std::string sdataset,
                        boost_array2d& array2D,
                        unsigned int nrow,
                        unsigned int ncol) {
  try
  {
    // Turn off the auto-printing
    H5::Exception::dontPrint();

    H5::Group group = file.openGroup(sgroup);

    H5::DataSet dataset = group.openDataSet(sdataset);

    H5::DataType type = dataset.getDataType();
//     /*
//       * Get the class of the datatype that is used by the dataset.
//       */
//     H5T_class_t type_class = dataset.getTypeClass();

    H5::DataSpace dataspace = dataset.getSpace();
    //Get the number of dimensions in the dataspace.
    int rank = dataspace.getSimpleExtentNdims();
    

    // Get the dimension size of each dimension in the dataspace
    hsize_t dims_out[1];    
    int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);

//     std::cerr << "rank " << rank << ", dimensions "
//               << static_cast<unsigned long>(dims_out[0]) << " x "
//               << static_cast<unsigned long>(dims_out[1]) << std::endl;

    // use boost multi array
    double* varray = new double[dims_out[0]*dims_out[1]];
    
    H5::DataSpace indataspace(rank, dims_out);
    
    dataset.read(&varray[0], H5::PredType::NATIVE_DOUBLE, indataspace, dataspace);
    
    for(unsigned int i = 0; i<dims_out[1]; ++i) {
//       std::cerr << std::endl;
      for(unsigned int j = 0; j<dims_out[0]; ++j) {
        array2D[i][j] = varray[i+nrow*j];
//         std::cerr << varray[i+nrow*j] << '\t';
      }
    }
    
    delete [] varray;
    // close dataspace and dataset
    indataspace.close();
    dataspace.close();    
    type.close();
    dataset.close();
    group.close();
   }  // end of try block
   // catch failure caused by the H5File operations
   catch( H5::FileIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataSet operations
   catch( H5::DataSetIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataSpace operations
   catch( H5::DataSpaceIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataSpace operations
   catch( H5::DataTypeIException error )
   {
      error.printError();
      return -1;
   }
   return 0;
}



template < class ArrType2D >
int write2d_hdf5(ArrType2D array2D, unsigned int nrow, unsigned int ncol,
                 std::string h5filename, std::string datasetname="data") {

//   H5std_string h5filename(filename.c_str());

//   const char* h5filename = filename.c_str();
//   std::cout << std::endl << "[ii] Try " << h5filename << std::endl;

  double* varray = new double[nrow*ncol];
    for(unsigned int i = 0; i<nrow; ++i) {
//       std::cout << std::endl;
      for(unsigned int j = 0; j<ncol; ++j) {
          varray[i+nrow*j] = array2D[i][j];
//           std::cout << varray[i+nrow*j] << '\t';
      }
  }
//   std::cout << std::endl << nrow << '\t' << ncol << std::endl;
  try
  {
    // Turn off the auto-printing
    H5::Exception::dontPrint();

    // Create new file
    H5::H5File file(h5filename, H5F_ACC_TRUNC);

    // dataset dimensions
    hsize_t dimsf[2];
    dimsf[0] = ncol;
    dimsf[1] = nrow;
    H5::DataSpace dataspace(2, dimsf);

    H5::DataType datatype(H5::PredType::NATIVE_DOUBLE);

    // FIXME datasetname is the same for all h5 files
    H5::DataSet dataset = file.createDataSet(datasetname, datatype, dataspace);

    //std::cout << std::endl << "[ii] Dataset filename: " << datasetname << std::endl;

    dataset.write(varray, H5::PredType::NATIVE_DOUBLE);

    delete [] varray;


    dataset.close();
    dataspace.close();
    file.close();
   }  // end of try block
   // catch failure caused by the H5File operations
   catch( H5::FileIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataSet operations
   catch( H5::DataSetIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataSpace operations
   catch( H5::DataSpaceIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataSpace operations
   catch( H5::DataTypeIException error )
   {
      error.printError();
      return -1;
   }
   return 0;
}

inline int open_hdf5(const std::string h5filename, H5::H5File& file,
                     std::string mode) {
  // ALERT some bytes lost. checked by valgrind = 1,360 bytes in 3 blocks
  try
  {
    // Turn off the auto-printing
    H5::Exception::dontPrint();

    if (mode == "create") {
      // Create new file / overwrite existing
      H5::H5File auxfile(h5filename, H5F_ACC_TRUNC);
      file = auxfile;
    }
    else {
      if (mode == "read") {
        H5::H5File auxfile(h5filename, H5F_ACC_RDONLY);
        file = auxfile;
      }
      else {
        if (mode == "write") {
          H5::H5File auxfile(h5filename, H5F_ACC_RDWR);
          file = auxfile;
        }
      }
    }
//     file = auxfile;

   }  // end of try block
   // catch failure caused by the H5File operations
   catch( H5::FileIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataSet operations
   catch( H5::DataSetIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataSpace operations
   catch( H5::DataSpaceIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataSpace operations
   catch( H5::DataTypeIException error )
   {
      error.printError();
      return -1;
   }
#ifdef VERBOSE
   std::cout << std::endl << "[ii] h5 file " << h5filename << " opened.\n";
#endif
   return 0;
}

inline int close_hdf5(H5::H5File& file) {

  try
  {
    // Turn off the auto-printing
    H5::Exception::dontPrint();

    // close dataspace and datafile
//     dataspace.close();
    file.close();

   }  // end of try block
   // catch failure caused by the H5File operations
   catch( H5::FileIException error )
   {
      error.printError();
      return -1;
   }
#ifdef VERBOSE
   std::cout << std::endl << "[ii] h5 file closed.\n";
#endif
   return 0;
}

template < class ArrType >
int livewrite_hdf5(H5::H5File& file, const ArrType array,
               const std::string dsname="data") {

  const unsigned int length = array.size();

  try
  {
    // Turn off the auto-printing
    H5::Exception::dontPrint();

    // WARNING could be defined in another function, overhead??
    // dataspace length
    // dataset dimensions
    hsize_t dimsf[1];
    dimsf[0] = length;
    H5::DataSpace dataspace(1, dimsf);

    // specify type double
    H5::DataType datatype(H5::PredType::NATIVE_DOUBLE);

    H5::DataSet dataset = file.createDataSet(dsname, datatype, dataspace);

    dataset.write(&array[0], H5::PredType::NATIVE_DOUBLE);

    // close dataspace and dataset
    dataspace.close();
    dataset.close();
   }  // end of try block
   // catch failure caused by the H5File operations
   catch( H5::FileIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataSet operations
   catch( H5::DataSetIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataSpace operations
   catch( H5::DataSpaceIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataSpace operations
   catch( H5::DataTypeIException error )
   {
      error.printError();
      return -1;
   }
   return 0;
}

template < class ArrType >
int liveread_hdf5(H5::H5File& file, ArrType &outarray,
               const std::string dsname="data") {

//   const unsigned int length = outarray.size();

  try
  {
    // Turn off the auto-printing
    H5::Exception::dontPrint();

    // WARNING could be defined in another function, overhead??
    // dataspace length
    // dataset dimensions
//     hsize_t dimsf[1];
//     dimsf[0] = length;
//     H5::DataSpace dataspace(1, dimsf);
//
//     // specify type double
//     H5::DataType datatype(H5::PredType::NATIVE_DOUBLE);
    H5::DataSet dataset = file.openDataSet( dsname );

    /*
      * Get the class of the datatype that is used by the dataset.
      */
    H5T_class_t type_class = dataset.getTypeClass();


    /*
      * Get dataspace of the dataset.
      */
    H5::DataSpace dataspace = dataset.getSpace();
    /*
      * Get the number of dimensions in the dataspace.
      */
    int rank = dataspace.getSimpleExtentNdims();
    /*
      * Get the dimension size of each dimension in the dataspace and
      * display them.
      */
    hsize_t dims_out[1];
    int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
//     std::cout << "rank " << rank << ", dimensions " <<
//         static_cast<unsigned long>(dims_out[0]) << " x ";/* <<
//         static_cast<unsigned long>(dims_out[1]) << std::endl;*/

    H5::DataSpace indataspace(1, dims_out);
    outarray.resize(static_cast<unsigned long>(dims_out[0]));
    dataset.read(&outarray[0], H5::PredType::NATIVE_DOUBLE, indataspace, dataspace);
//     for (unsigned int i = 0; i<outarray.size(); ++i)
//     {
//       std::cout << std::endl << outarray[i];
//     }
    // close dataspace and dataset
    dataspace.close();
    dataset.close();
   }  // end of try block
   // catch failure caused by the H5File operations
   catch( H5::FileIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataSet operations
   catch( H5::DataSetIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataSpace operations
   catch( H5::DataSpaceIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataSpace operations
   catch( H5::DataTypeIException error )
   {
      error.printError();
      return -1;
   }
   return 0;
}

// create unsigned attribute
inline
int create_attrib_hdf5(H5::H5File& file, const std::string sgroup,
                      const std::string sattr, const unsigned int &value) {
  try
  {
    // Turn off the auto-printing
    H5::Exception::dontPrint();

    H5::Group group;
    // attempt to open group
    try {
      H5::Exception::dontPrint();
      group = file.openGroup(sgroup);
    }
    catch( H5::Exception gerror ) {
//       gerror.printError();
      // if group does not exist, create it
      try {
        H5::Exception::dontPrint();
        group = file.createGroup(sgroup);
      }
      catch( H5::Exception gerror2 ) {
        gerror2.printError();
        std::cout << std::endl << "[ee] h5 I/O error in group "
                  << sgroup << ". Terminate.\n";
        std::terminate();
      }
    }

    H5::Attribute attr;
    // attempt to open group
    try {
      H5::Exception::dontPrint();
      attr = group.openAttribute(sattr);
    }
    catch( H5::Exception aerror )
    {
//       aerror.printError();
      try {
        H5::DataSpace dspace(H5S_SCALAR);
        H5::DataType type(H5::PredType::NATIVE_UINT);
        attr = group.createAttribute(sattr, type, dspace);
        type.close();
        dspace.close();
      }
      catch( H5::Exception aerror2 )
      {
        aerror2.printError();
        std::cout << std::endl << "[ee] h5 I/O error in attribute "
                  << sattr << ". Terminate.\n";
        std::terminate();
      }
    }
    H5::DataType type = attr.getDataType();
//
    void *buf = nullptr;
//
    std::memcpy(&buf, &value, sizeof(value));
//
    attr.write(type, &buf);
//
    type.close();
//
    attr.close();

    group.close();
   }  // end of try block
   // catch failure caused by the H5File operations
   catch( H5::FileIException error )
   {
      error.printError();
      std::cout << std::endl << "[ee] h5 I/O error in file. Terminate.\n";
      std::terminate();
      return -1;
   }
   // catch failure caused by the Group operations
   catch( H5::GroupIException error )
   {
      error.printError();
      std::cout << std::endl << "[ee] h5 I/O error in group "
                << sgroup << ". Terminate.\n";
      std::terminate();
      return -1;
   }
   // catch failure caused by the Attribute operations
   catch( H5::AttributeIException error )
   {
      error.printError();
      std::cout << std::endl << "[ee] h5 I/O error in attribute "
                << sattr << ". Terminate.\n";
      std::terminate();
      return -1;
   }
   // catch failure caused by the DataType operations
   catch( H5::DataTypeIException error )
   {
      error.printError();
      std::cout << std::endl << "[ee] h5 I/O error in group DataType. Terminate.\n";
      std::terminate();
      return -1;
   }
   return 0;
}


// TODO use template see below create_dsattrib_hdf5
// create double attribute
inline
int create_attrib_hdf5(H5::H5File& file, const std::string sgroup,
                      const std::string sattr, const double &value) {
  try
  {
    // Turn off the auto-printing
    H5::Exception::dontPrint();

    H5::Group group;
    // attempt to open group
    try {
      H5::Exception::dontPrint();
      group = file.openGroup(sgroup);
    }
    catch( H5::Exception gerror ) {
//       gerror.printError();
      // if group does not exist, create it
      try {
        H5::Exception::dontPrint();
        group = file.createGroup(sgroup);
      }
      catch( H5::Exception gerror2 ) {
        gerror2.printError();
        std::cout << std::endl << "[ee] h5 I/O error in group "
                  << sgroup << ". Terminate.\n";
        std::terminate();
      }
    }

    H5::Attribute attr;
    // attempt to open group
    try {
      H5::Exception::dontPrint();
      attr = group.openAttribute(sattr);
    }
    catch( H5::Exception aerror )
    {
//       aerror.printError();
      try {
        H5::DataSpace dspace(H5S_SCALAR);
        H5::DataType type(H5::PredType::NATIVE_DOUBLE);
        attr = group.createAttribute(sattr, type, dspace);
        type.close();
        dspace.close();
      }
      catch( H5::Exception aerror2 )
      {
        aerror2.printError();
        std::cout << std::endl << "[ee] h5 I/O error in attribute "
                  << sattr << ". Terminate.\n";
        std::terminate();
      }
    }
    H5::DataType type = attr.getDataType();
//
    void *buf = nullptr;
//
    std::memcpy(&buf, &value, sizeof(value));
//
    attr.write(type, &buf);
//
    type.close();
//
    attr.close();

    group.close();
   }  // end of try block
   // catch failure caused by the H5File operations
   catch( H5::FileIException error )
   {
      error.printError();
      std::cout << std::endl << "[ee] h5 I/O error in file. Terminate.\n";
      std::terminate();
      return -1;
   }
   // catch failure caused by the Group operations
   catch( H5::GroupIException error )
   {
      error.printError();
      std::cout << std::endl << "[ee] h5 I/O error in group "
                << sgroup << ". Terminate.\n";
      std::terminate();
      return -1;
   }
   // catch failure caused by the Attribute operations
   catch( H5::AttributeIException error )
   {
      error.printError();
      std::cout << std::endl << "[ee] h5 I/O error in attribute "
                << sattr << ". Terminate.\n";
      std::terminate();
      return -1;
   }
   // catch failure caused by the DataType operations
   catch( H5::DataTypeIException error )
   {
      error.printError();
      std::cout << std::endl << "[ee] h5 I/O error in group DataType. Terminate.\n";
      std::terminate();
      return -1;
   }
   return 0;
}

// overloading types
inline
void set_datatype(H5::DataType& type, double value){
  type = H5::DataType(H5::PredType::NATIVE_DOUBLE);
} 

inline
void set_datatype(H5::DataType& type, unsigned int value){
  type = H5::DataType(H5::PredType::NATIVE_UINT);
} 

inline
void set_datatype(H5::DataType& type, int value){
  type = H5::DataType(H5::PredType::NATIVE_INT);
} 

inline
void set_datatype(H5::DataType& type, long int value){
  type = H5::DataType(H5::PredType::NATIVE_LONG);
} 

template < class Type >
int create_dsattrib_hdf5(H5::H5File& file, const std::string sgroup,
                         const std::string dsname, const std::string sattr,
                         const Type &value) {
  try
  {
    // Turn off the auto-printing
    H5::Exception::dontPrint();

    H5::Group group;
    H5::DataSet dataset;
    // attempt to open group

    group = file.openGroup(sgroup);

    dataset = group.openDataSet(dsname);

    H5::Attribute attr;
 
    H5::DataSpace dspace(H5S_SCALAR);
//     H5::DataType type(H5::PredType::NATIVE_DOUBLE);
    H5::DataType type;
    set_datatype(type, value);

    attr = dataset.createAttribute(sattr, type, dspace);       

//
    void *buf = nullptr;
//
    std::memcpy(&buf, &value, sizeof(value));
//
    attr.write(type, &buf);
//
    type.close();
//
    dspace.close();
  
    attr.close();

    dataset.close();
    
    group.close();
   }  // end of try block
   // catch failure caused by the H5File operations
   catch( H5::FileIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the Group operations
   catch( H5::GroupIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataSet operations
   catch( H5::DataSetIException error )
   {
      error.printError();
      return H5_DSET_ERR;
   }
   // catch failure caused by the DataSpace operations
   catch( H5::DataSpaceIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataSpace operations
   catch( H5::DataTypeIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataSet operations
   catch( H5::AttributeIException error )
   {
      error.printError();
      return -1;
   }
   return 0;
}

template < class Type >
int read_dsattrib_hdf5(H5::H5File& file, const std::string sgroup,
                       const std::string dsname, const std::string sattr,
                       Type &value) {
  try
  {
    // Turn off the auto-printing
    H5::Exception::dontPrint();

    H5::Group group = file.openGroup(sgroup);
    
    H5::DataSet dataset  = group.openDataSet(dsname);
    // attempt to open group

    H5::Attribute attr = dataset.openAttribute(sattr);

    H5::DataType type = attr.getDataType();
    
    void *buf = nullptr;
    attr.read(type, &buf);

//     assert(sizeof(value) == sizeof(buf));
    std::memcpy(&value, &buf, sizeof(value));

    type.close();

    attr.close();
   
    group.close();    
   }  // end of try block
   // catch failure caused by the H5File operations
   catch( H5::FileIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the Group operations
   catch( H5::GroupIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataSet operations
   catch( H5::DataSetIException error )
   {
      error.printError();
      return H5_DSET_ERR;
   }
   // catch failure caused by the DataSpace operations
   catch( H5::DataSpaceIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataSpace operations
   catch( H5::DataTypeIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataSet operations
   catch( H5::AttributeIException error )
   {
      error.printError();
      return -1;
   }
   return 0;
}



template < class Type >
int write_attrib_hdf5(H5::H5File& file, const std::string sgroup,
                     const std::string sattr, Type &value) {
  try
  {
    // Turn off the auto-printing
    H5::Exception::dontPrint();

    H5::Group group = file.openGroup(sgroup);

    H5::Attribute attr = group.openAttribute(sattr);

    H5::DataType type = attr.getDataType();

    void *buf = nullptr;

    std::memcpy(&buf, &value, sizeof(value));

    attr.write(type, &buf);

    type.close();

    attr.close();

    group.close();
   }  // end of try block
   // catch failure caused by the H5File operations
   catch( H5::FileIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the Group operations
   catch( H5::GroupIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the Attribute operations
   catch( H5::AttributeIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataType operations
   catch( H5::DataTypeIException error )
   {
      error.printError();
      return -1;
   }
   return 0;
}

template < class Type >
int read_attrib_hdf5(H5::H5File& file, const std::string sgroup,
                     const std::string sattr, Type &value) {
  try
  {
    // Turn off the auto-printing
    H5::Exception::dontPrint();

    H5::Group group = file.openGroup(sgroup);

    H5::Attribute attr = group.openAttribute(sattr);

    H5::DataType type = attr.getDataType();

    void *buf = nullptr;
    attr.read(type, &buf);

//     assert(sizeof(value) == sizeof(buf));
    std::memcpy(&value, &buf, sizeof(value));

    type.close();

    attr.close();

    group.close();
   }  // end of try block
   // catch failure caused by the H5File operations
   catch( H5::FileIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the Group operations
   catch( H5::GroupIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the Attribute operations
   catch( H5::AttributeIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataType operations
   catch( H5::DataTypeIException error )
   {
      error.printError();
      return -1;
   }
   return 0;
}

// Read string Attribute
inline
int read_str_attrib_hdf5(H5::H5File& file, const std::string sgroup,
                        const std::string sattr, std::string& stdstr) {
  try
  {
    // Turn off the auto-printing
    H5::Exception::dontPrint();

    H5::Group group = file.openGroup(sgroup);

    H5::Attribute attr = group.openAttribute(sattr);

    H5::DataType type = attr.getDataType();

//     H5std_string h5str;

    attr.read(type, stdstr);

    type.close();

    attr.close();

    group.close();
   }  // end of try block
   // catch failure caused by the H5File operations
   catch( H5::FileIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the Group operations
   catch( H5::GroupIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the Attribute operations
   catch( H5::AttributeIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataType operations
   catch( H5::DataTypeIException error )
   {
      error.printError();
      return -1;
   }
   return 0;
}

// WARNING reading doubles
// Read dataset in group
template < class ArrType >
int read_dataset_hdf5(H5::H5File& file, ArrType &outarray,
                      const std::string sgroup,
                      const std::string dsname) {

  try
  {
    // Turn off the auto-printing
    H5::Exception::dontPrint();

    H5::Group group = file.openGroup(sgroup);

    H5::DataSet dataset = group.openDataSet( dsname );

    /*
     * Get the class of the datatype that is used by the dataset.
     */
    H5T_class_t type_class = dataset.getTypeClass();
    /*
     * Get dataspace of the dataset.
     */
    H5::DataSpace dataspace = dataset.getSpace();
    /*
     * Get the number of dimensions in the dataspace.
     */
    int rank = dataspace.getSimpleExtentNdims();
    /*
     * Get the dimension size of each dimension in the dataspace.
     */
    hsize_t dims_out[1];
    int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
//     std::cout << "rank " << rank << ", dimensions " <<
//         static_cast<unsigned long>(dims_out[0]) << " x ";/* <<
//         static_cast<unsigned long>(dims_out[1]) << std::endl;*/

    H5::DataSpace indataspace(1, dims_out);
    outarray.resize(static_cast<unsigned long>(dims_out[0]));
    dataset.read(&outarray[0], H5::PredType::NATIVE_DOUBLE, indataspace, dataspace);
//     for (unsigned int i = 0; i<outarray.size(); ++i)
//     {
//       std::cout << std::endl << outarray[i];
//     }
    // close dataspace and dataset
    indataspace.close();
    
    dataspace.close();

    dataset.close();

    group.close();
   }  // end of try block
   // catch failure caused by the H5File operations
   catch( H5::FileIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the Group operations
   catch( H5::GroupIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataSet operations
   catch( H5::DataSetIException error )
   {
      error.printError();
      return H5_DSET_ERR;
   }
   // catch failure caused by the DataSpace operations
   catch( H5::DataSpaceIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataSpace operations
   catch( H5::DataTypeIException error )
   {
      error.printError();
      return -1;
   }
   return 0;
}

// write dataset in group
template < class ArrType >
int write_dataset_hdf5(H5::H5File& file, const ArrType &outarray,
                      const std::string sgroup,
                      const std::string dsname) {

  try
  {
    // Turn off the auto-printing
    H5::Exception::dontPrint();

    H5::Group group = file.openGroup(sgroup);

    H5::DataSet dataset = group.openDataSet( dsname );

    /*
     * Get the class of the datatype that is used by the dataset.
     */
    H5T_class_t type_class = dataset.getTypeClass();
    /*
     * Get dataspace of the dataset.
     */
    H5::DataSpace dataspace = dataset.getSpace();
    /*
     * Get the number of dimensions in the dataspace.
     */
    int rank = dataspace.getSimpleExtentNdims();
    /*
     * Get the dimension size of each dimension in the dataspace.
     */
    hsize_t dims_out[1];
    int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
//     std::cout << "rank " << rank << ", dimensions " <<
//         static_cast<unsigned long>(dims_out[0]) << " x ";/* <<
//         static_cast<unsigned long>(dims_out[1]) << std::endl;*/

    H5::DataSpace indataspace(1, dims_out);
//     outarray.resize(static_cast<unsigned long>(dims_out[0]));
//     dataset.read(&outarray[0], H5::PredType::NATIVE_DOUBLE, indataspace, dataspace);
    dataset.write(&outarray[0], H5::PredType::NATIVE_DOUBLE);

    // close dataspace and dataset
    dataspace.close();

    dataset.close();

    group.close();
   }  // end of try block
   // catch failure caused by the H5File operations
   catch( H5::FileIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the Group operations
   catch( H5::GroupIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataSet operations
   catch( H5::DataSetIException error )
   {
      error.printError();
      return H5_DSET_ERR;
   }
   // catch failure caused by the DataSpace operations
   catch( H5::DataSpaceIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataSpace operations
   catch( H5::DataTypeIException error )
   {
      error.printError();
      return -1;
   }
   return 0;
}

template < class ArrType >
int create_dataset_hdf5(H5::H5File& file, const ArrType &outarray,
                        const std::string sgroup,
                        const std::string dsname) {

  try
  {
    // Turn off the auto-printing
    H5::Exception::dontPrint();

    H5::Group group;
    // attempt to open group
    try {
      H5::Exception::dontPrint();
      group = file.openGroup(sgroup);
    }
    catch( H5::Exception error_opengroup ) {
//       error_opengroup.printError();
      // if group does not exist, create it
      try {
        H5::Exception::dontPrint();
        group = file.createGroup(sgroup);
      }
      catch( H5::Exception error_creategroup ) {
        error_creategroup.printError();
        std::cout << std::endl << "[ee] h5 I/O error in group "
                  << sgroup << ". Terminate.\n";
        std::terminate();
      }
    }

    H5::DataSet dataset;
    // attempt to open dataset
    try {
      H5::Exception::dontPrint();
      dataset = group.openDataSet(dsname);
    }
    catch( H5::Exception error_opendset ) {
//       error_opendset.printError();
      // if dataset not exist, create it
      try {
        H5::Exception::dontPrint();
        // specify type double
        H5::DataType datatype(H5::PredType::NATIVE_DOUBLE);

        hsize_t dimsf[1];
        dimsf[0] = outarray.size();
        H5::DataSpace dataspace(1, dimsf);
        dataset = group.createDataSet(dsname, datatype, dataspace);
        datatype.close();
        dataspace.close();
      }
      catch( H5::Exception error_createdset )
      {
        error_createdset.printError();
        std::cout << std::endl << "[ee] h5 I/O error in dataset "
                  << dsname << ". Terminate.\n";
        std::terminate();
      }
    }
    H5::DataType datatype(H5::PredType::NATIVE_DOUBLE);

    dataset.write(&outarray[0], datatype);

    // close dataspace and dataset
    datatype.close();
    dataset.close();
//     group.close();

   }  // end of try block
   // catch failure caused by the H5File operations
   catch( H5::FileIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the Group operations
   catch( H5::GroupIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataSet operations
   catch( H5::DataSetIException error )
   {
      error.printError();
      return H5_DSET_ERR;
   }
   // catch failure caused by the DataSpace operations
   catch( H5::DataSpaceIException error )
   {
      error.printError();
      return -1;
   }
   // catch failure caused by the DataSpace operations
   catch( H5::DataTypeIException error )
   {
      error.printError();
      return -1;
   }
   return 0;
}

#endif/* H5PLASMA_H*/
