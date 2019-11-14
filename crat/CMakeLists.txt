
cmake_minimum_required(VERSION 3.0)

project(crat)

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

# defines variable SRCS to include cpp files
file(GLOB SRCS *.cpp)
file(GLOB INCSRCS RELATIVE "../include/" "*.cpp")

include_directories("../include/")
include_directories(SYSTEM "/usr/include/eigen3")

add_definitions(-DH5_USE_18_API)
add_definitions(-DBOOST_LOG_DYN_LINK)
add_definitions(-DNDEBUG)
add_definitions(-DBOOST_UBLAS_DNDEBUG)
add_definitions(-DUSING_EIGEN)
add_definitions(-DEIGEN_DONT_PARALLELIZE)

add_library(GridModel STATIC ../include/GridModel.cpp)
add_library(CRate STATIC ../include/CRate.cpp)

# Configure ndust version
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/../include/ndust_version.cpp.in" "${CMAKE_CURRENT_SOURCE_DIR}/../include/ndust_version.cpp" @ONLY)
add_library(ndust_version STATIC ../include/ndust_version.cpp)

# GSL
set(LIBS "-lgsl -lgslcblas -lm -llapack -lumfpack -lboost_log -lboost_thread -lboost_system -lboost_log_setup")

# HDF5
set (LIBS "${LIBS} -I/usr/include/hdf5/ -lhdf5_cpp -lhdf5")
set (LIBS "${LIBS} -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_BSD_SOURCE -D_FORTIFY_SOURCE=2 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Werror=format-security -lpthread -lz -ldl -lm")

# include_directories(${INC_DIR})

# -frecord-gcc-switches
# records flags passed to gcc
# can be retrieved by $readelf -p .GCC.command.line crat
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=gnu++14 -march=znver1 -Ofast -ftree-vectorize -fopt-info-vec-missed=vec.miss -frecord-gcc-switches -fgraphite-identity -ftree-loop-distribution -floop-nest-optimize -flto")

set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} -Ofast -frecord-gcc-switches -march=znver1")

set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -Ofast -frecord-gcc-switches -march=znver1")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=gnu++14 -march=znver1 -O3 -ftree-vectorize -fopt-info-vec-missed=vec.miss -frecord-gcc-switches")
#
# RELEASE ADD -DNDEBUG for boost performance
#
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DNDEBUG -fopenmp -DBOOST_LOG_DYN_LINK")
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

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -frecord-gcc-switches -march=znver1")

add_executable(crat "main.cpp" ${INCSRCS})

add_executable(crat-writer "crat-writer.cpp" ${INCSRCS})

add_executable(crat-reader "crat-reader.cpp" ${INCSRCS})

add_executable(crat-merge "crat-merge.cpp" ${INCSRCS})

add_executable(potentials "potentials.cpp" ${INCSRCS})

add_executable(potential_b1 "potentials_b1.cpp" ${INCSRCS})

target_link_libraries(CRate
                      GridModel
                      ${LIBS})

target_link_libraries(crat
                      CRate
		      ndust_version
                      ${LIBS})

target_link_libraries(crat-writer
                      CRate
		      ndust_version
                      ${LIBS})

target_link_libraries(crat-reader
                      CRate
		      ndust_version
                      ${LIBS})

target_link_libraries(crat-merge
                      CRate
		      ndust_version
                      ${LIBS})

target_link_libraries(potentials ${LIBS})

target_link_libraries(potential_b1 ${LIBS})

enable_testing()

install(TARGETS crat RUNTIME DESTINATION bin)

