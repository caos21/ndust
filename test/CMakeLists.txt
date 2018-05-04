cmake_minimum_required(VERSION 3.0)

set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 

add_definitions(-DBOOST_ALL_NO_LIB)

include_directories(SYSTEM "/usr/include/eigen3")

add_definitions(-DBOOST_LOG_DYN_LINK)
#add_definitions(-DNDEBUG)
#add_definitions(-DBOOST_UBLAS_DNDEBUG)
#add_definitions(-DUSING_EIGEN)
#add_definitions(-DEIGEN_DONT_PARALLELIZE)


#set(LIBS "-lgsl -lgslcblas -lm -llapack -lumfpack -lboost_log -lboost_thread -lboost_system -lboost_log_setup")

set(LIBS "-llapack")

find_package(Boost COMPONENTS log thread system atomic chrono date_time unit_test_framework REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

#message("${Boost_LIBRARIES}")
#link_directories("${BOOST}/lib/")

find_package(OpenMP)
if(OPENMP_FOUND)
 message("** OpenMP found")
 set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
 set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} ")
 set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS} ")
 add_definitions(-DOPENMP)
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O1 -frecord-gcc-switches -Wall -pedantic -march=znver1 -g")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LIBS}")

add_executable(eint-test eint-test.cpp)

target_link_libraries(eint-test ${Boost_LIBRARIES} ${LIBS})

add_executable(enhancement-test enhancement-test.cpp)

target_link_libraries(enhancement-test ${Boost_LIBRARIES} ${LIBS})
