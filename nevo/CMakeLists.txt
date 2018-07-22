
cmake_minimum_required(VERSION 3.0)

project(nevo)

# # defines variable SRCS to include cpp files
# file(GLOB SRCS *.cpp)
# file(GLOB INCSRCS RELATIVE "../include/" "*.cpp")

set(CVODE_LIB sundials_cvode)
set(NVECS_LIB sundials_nvecserial)
set(EXTRA_LINK_LIBS m)

set(SUNDIALS_LIBS ${CVODE_LIB} ${NVECS_LIB} ${EXTRA_LINK_LIBS})

include_directories("../include/")

link_directories("../build/")

add_library(GridModel STATIC ../include/GridModel.cpp)
add_library(CRate STATIC ../include/CRate.cpp)
add_library(PlasmaModel STATIC ../include/PlasmaModel.cpp)
add_library(NanoModel STATIC ../include/NanoModel.cpp)
add_library(NEvo STATIC ../NEvo.cpp)

# GSL
set(LIBS "-lgsl -lgslcblas -lm -llapack -lumfpack -lboost_log -lboost_thread -lboost_system -lboost_log_setup -lhdf5_cpp -lhdf5 -lpthread -lz -ldl -lm")

set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} -ftree-vectorizer-verbose=3 -std=gnu++14 -march=znver1 -O3 -ftree-vectorize -fstack-protector --param=ssp-buffer-size=4 -Wformat -Werror=format-security")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ftree-vectorizer-verbose=3 -std=gnu++14 -march=znver1 -O3 -ftree-vectorize -fopt-info-vec-missed=vec.miss -fstack-protector --param=ssp-buffer-size=4 -Wformat -Werror=format-security")

# RELEASE ADD -DNDEBUG for boost performance
add_definitions(-DBOOST_LOG_DYN_LINK -DNDEBUG -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_DEFAULT_SOURCE)


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
                      ${LIBS})

target_link_libraries(nevo
                      NEvo
                      CRate
                      PlasmaModel
                      NanoModel
                      ${LIBS}
                      ${SUNDIALS_LIBS})

enable_testing()

install(TARGETS nevo RUNTIME DESTINATION bin)


