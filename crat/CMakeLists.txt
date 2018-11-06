
cmake_minimum_required(VERSION 2.8)

project(crat)

set(CMAKE_MODULE_RELPATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")

message("${CMAKE_MODULE_RELPATH}")
#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_MODULE_RELPATH}")

include(GetGitRevisionDescription)

get_git_head_revision(NDUST_BRANCH NDUST_HASH)
message(${NDUST_HASH})
message(${NDUST_BRANCH})

git_local_changes(NDUST_STATE)
message(${NDUST_STATE})

git_describe(NDUST_VERSION "--long")
message(${NDUST_VERSION})


include_directories(SYSTEM "/usr/include/eigen3")

add_definitions(-DBOOST_LOG_DYN_LINK)
add_definitions(-DNDEBUG)
add_definitions(-DBOOST_UBLAS_DNDEBUG)
add_definitions(-DUSING_EIGEN)
add_definitions(-DEIGEN_DONT_PARALLELIZE)


set (CMAKE_CXX_FLAGS "-std=gnu++14 ${CMAKE_CXX_FLAGS}")

find_package(Boost COMPONENTS thread system atomic chrono date_time REQUIRED)

# defines variable SRCS to include cpp files
file(GLOB SRCS *.cpp)
file(GLOB INCSRCS RELATIVE "../include/" "*.cpp")

# project includes
include_directories("../include/")
# include local boost-numeric-bindings
include_directories("~/local/include/boost-numeric-bindings/")
# include  boost
#include_directories("~/lib/boost/include/")
# include lapack blas cblas
#include_directories("/home/bsantos/local/include")
# include GSL
include_directories("/cvmfs/soft.computecanada.ca/easybuild/software/2017/sse3/Compiler/intel2016.4/gsl/2.3/include/")
# include HDF5
include_directories("/cvmfs/soft.computecanada.ca/easybuild/software/2017/sse3/Compiler/intel2016.4/hdf5/1.8.18/include/")

# local boost
#link_directories("~/lib/boost/lib/")
# lapack blas cblas
#link_directories("/home/bsantos/local/lib/")
# gsl
link_directories("/cvmfs/soft.computecanada.ca/easybuild/software/2017/sse3/Compiler/intel2016.4/gsl/2.3/lib/")
# HDF5
link_directories("/cvmfs/soft.computecanada.ca/easybuild/software/2017/sse3/Compiler/intel2016.4/hdf5/1.8.18/lib/")


add_library(GridModel STATIC ../include/GridModel.cpp)
add_library(CRate STATIC ../include/CRate.cpp)


# include_directories(${INC_DIR})


set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} -Ofast -frecord-gcc-switches")

set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -Ofast -frecord-gcc-switches")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=gnu++14 -mtune=native -O2 -frecord-gcc-switches")

#
# RELEASE ADD -DNDEBUG for boost performance
#
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DNDEBUG -fopenmp")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fopenmp")

list(APPEND ${CMAKE_CXX_COMPILE_FLAGS} ${LIBS})

# WARNING not working
find_package(OpenMP)
if(OPENMP_FOUND)
  message("** OpenMP found")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} ")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS} ")
  add_definitions(-DOPENMP)
endif()


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -frecord-gcc-switches")

# Configure ndust version
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/../include/ndust_version.cpp.in" "${CMAKE_CURRENT_SOURCE_DIR}/../include/ndust_version.cpp" @ONLY)

add_library(ndust_version STATIC ../include/ndust_version.cpp)

add_executable(crat "main.cpp" ${INCSRCS})

add_executable(crat-writer "crat-writer.cpp" ${INCSRCS})

add_executable(crat-reader "crat-reader.cpp" ${INCSRCS})

add_executable(crat-merge "crat-merge.cpp" ${INCSRCS})

# target libraries
#list(APPEND TLIBS "boost_log hdf5 hdf5_cpp boost_system boost_thread boost_atomic boost_chrono boost_regex boost_date_time boost_filesystem mkl")

target_link_libraries(CRate
  GridModel
  boost_log
  ndust_version
  hdf5
  hdf5_cpp
  boost_system
  boost_thread
  boost_atomic
  boost_chrono
  boost_regex
  boost_date_time
  boost_filesystem
  mkl)

target_link_libraries(crat
  CRate
  boost_log
  hdf5
  hdf5_cpp
  boost_system
  boost_thread
  boost_atomic
  boost_chrono
  boost_regex
  boost_date_time
  boost_filesystem
  mkl)

target_link_libraries(crat-writer
  CRate
  boost_log
  hdf5
  hdf5_cpp
  boost_system
  boost_thread
  boost_atomic
  boost_chrono
  boost_regex
  boost_date_time
  boost_filesystem
  mkl)

target_link_libraries(crat-reader
  CRate
  boost_log
  hdf5
  hdf5_cpp
  boost_system
  boost_thread
  boost_atomic
  boost_chrono
  boost_regex
  boost_date_time
  boost_filesystem
  mkl)

target_link_libraries(crat-merge
  CRate
  boost_log
  hdf5
  hdf5_cpp
  boost_system
  boost_thread
  boost_atomic
  boost_chrono
  boost_regex
  boost_date_time
  boost_filesystem
  mkl)
  

enable_testing()

install(TARGETS crat RUNTIME DESTINATION bin)

