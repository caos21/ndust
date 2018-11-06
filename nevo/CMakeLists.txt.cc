
cmake_minimum_required(VERSION 2.8)

project(nevo)

set(CMAKE_MODULE_RELPATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")

message("oo cmake module relative path : ${CMAKE_MODULE_RELPATH}")
#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_MODULE_RELPATH}")

include(GetGitRevisionDescription)

get_git_head_revision(NDUST_BRANCH NDUST_HASH)
message("oo Hash : ${NDUST_HASH}")
message("oo Origin : ${NDUST_BRANCH}")

git_local_changes(NDUST_STATE)
message("oo State : ${NDUST_STATE}")

git_describe(NDUST_VERSION "--long")
message("oo Version ${NDUST_VERSION}")


# ============= LSODA
get_filename_component(PARENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)
set(LSODA_DIR "${PARENT_DIR}/modules/liblsoda/")
include(ExternalProject)
ExternalProject_Add(liblsoda
    SOURCE_DIR ${LSODA_DIR}
    CONFIGURE_COMMAND ""
    STEP_TARGETS build
    BUILD_COMMAND make
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND "")
include_directories("${LSODA_DIR}/src/")
link_directories("${LSODA_DIR}/src/")
# ============

add_definitions(-DBOOST_LOG_DYN_LINK)
add_definitions(-DNDEBUG)
add_definitions(-DBOOST_UBLAS_DNDEBUG)
add_definitions(-DEIGEN_DONT_PARALLELIZE)

# Specify path to SUNDIALS header files
SET(SUNDIALS_INC_DIR
  /home/bsantos/local/sundials3/include
  CACHE STRING
  "Location of SUNDIALS header files")

# Add path to SUNDIALS header files
INCLUDE_DIRECTORIES(${SUNDIALS_INC_DIR})

# Set search path for SUNDIALS libraries 
SET(SUNDIALS_LIB_DIR /home/bsantos/local/sundials3/lib)

# Find the SUNDIALS solvers library
FIND_LIBRARY(SUNDIALS_SOLVER_LIB
  sundials_cvode ${SUNDIALS_LIB_DIR}
  DOC "CVODE library")

# Find the NVECTOR library
FIND_LIBRARY(SUNDIALS_NVEC_LIB
  sundials_nvecserial ${SUNDIALS_LIB_DIR}
  DOC "NVECTOR library")


# Set additional libraries
SET(SUNDIALS_EXTRA_LIB  -lm /cvmfs/soft.computecanada.ca/nix/var/nix/profiles/16.09/lib64/librt.so CACHE STRING "Additional libraries")

# List of Sundials libraries shared across all examples
SET(SUNDIALS_LIBS ${SUNDIALS_SOLVER_LIB} ${SUNDIALS_NVEC_LIB} ${SUNDIALS_EXTRA_LIB})

set (CMAKE_CXX_FLAGS "-std=gnu++14 ${CMAKE_CXX_FLAGS}")

include_directories("../include/")

#link_directories("../build/")

add_library(GridModel STATIC ../include/GridModel.cpp)
add_library(CRate STATIC ../include/CRate.cpp)
add_library(PlasmaModel STATIC ../include/PlasmaModel.cpp)
add_library(NanoModel STATIC ../include/NanoModel.cpp)
add_library(NEvo STATIC ../NEvo.cpp)

# Configure ndust version
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/../include/ndust_version.cpp.in" "${CMAKE_CURRENT_SOURCE_DIR}/../include/ndust_version.cpp" @ONLY)
add_library(ndust_version STATIC ../include/ndust_version.cpp)

# include local boost-numeric-bindings
include_directories("~/local/include/boost-numeric-bindings/")
# include GSL
include_directories("/cvmfs/soft.computecanada.ca/easybuild/software/2017/sse3/Compiler/intel2016.4/gsl/2.3/include/")
# include HDF5
include_directories("/cvmfs/soft.computecanada.ca/easybuild/software/2017/sse3/Compiler/intel2016.4/hdf5/1.8.18/include/")

#link_directories("/home/bsantos/local/lib/")
# gsl
link_directories("/cvmfs/soft.computecanada.ca/easybuild/software/2017/sse3/Compiler/intel2016.4/gsl/2.3/lib/")
# HDF5
link_directories("/cvmfs/soft.computecanada.ca/easybuild/software/2017/sse3/Compiler/intel2016.4/hdf5/1.8.18/lib/")


set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} -Ofast -frecord-gcc-switches")

set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -Ofast -frecord-gcc-switches")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=gnu++14 -Og")

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



find_package(HDF5)
if(HDF5_VERSION VERSION_LESS "1.10.0")
  # use old api 1.8
  message("HDF5 version 1.10 < ${HDF5_VERSION}")
else()
  # use new api
  message("HDF5 version ${HDF5_VERSION} > 1.10.0")
  add_definitions("-DH5API110")
endif()


add_executable(nevo main.cpp)

target_link_libraries(NEvo
                      GridModel
                      NanoModel
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
		      lsoda
		      mkl
		      ${LIBS}
		      ${SUNDIALS_LIBS})

target_link_libraries(nevo
                      NEvo
                      CRate
                      PlasmaModel
                      NanoModel
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
		      lsoda
		      mkl
                      ${LIBS}
                      ${SUNDIALS_LIB})

enable_testing()

install(TARGETS nevo RUNTIME DESTINATION bin)



