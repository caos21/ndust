
cmake_minimum_required(VERSION 3.0)

project(crat)

# defines variable SRCS to include cpp files
file(GLOB SRCS *.cpp)
file(GLOB INCSRCS RELATIVE "../include/" "*.cpp")

include_directories("../include/")

add_library(GridModel STATIC ../include/GridModel.cpp)
add_library(CRate STATIC ../include/CRate.cpp)

# GSL
set(LIBS "-lgsl -lgslcblas -lm -llapack -lumfpack -lboost_log -lboost_thread -lboost_system -lboost_log_setup")

# HDF5
set (LIBS "${LIBS} -I/usr/include/hdf5/ -lhdf5_cpp -lhdf5")
set (LIBS "${LIBS} -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_BSD_SOURCE -D_FORTIFY_SOURCE=2 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Werror=format-security -lpthread -lz -ldl -lm")

# include_directories(${INC_DIR})

set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} -O2 -ftree-vectorizer-verbose=3 ")

set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=gnu++11 -mtune=athlon-fx -O2 ")
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

add_executable(crat ${SRCS} ${INCSRCS})

target_link_libraries(CRate
                      GridModel
                      ${LIBS})

target_link_libraries(crat                      
                      CRate
                      ${LIBS})

enable_testing()

install(TARGETS crat RUNTIME DESTINATION bin)

