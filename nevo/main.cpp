#include <iostream>
#include <string>

// OpenMP
#include <omp.h>

// io functions
#include "io.h"

// logger
#include "log.h"

// class CRate
// #include "CRate.h"

#include "NEvo.h"

int main(int argc, char **argv) {
  // Prints header
  print_header(std::cout);

  // Print help
  show_help(argc, argv);

  // Checks if option -g is called
  char *input_opt;
  input_opt = get_cmd_option(argv, argv + argc, "-g");
  std::string grid_filename = "";
  if(input_opt) {
    // sets new grid filename
    grid_filename.append(input_opt);
  }
 
  // Checks if option -p is called
  input_opt = get_cmd_option(argv, argv + argc, "-p");
  std::string plasma_filename = "";
  if(input_opt) {
    // sets new plasma filename
    plasma_filename.append(input_opt);
  }

  // Checks if option -n is called
  input_opt = get_cmd_option(argv, argv + argc, "-n");
  std::string nano_filename = "";
  if(input_opt) {
    // sets new plasma filename
    nano_filename.append(input_opt);
  }
  
  
  // Check environment variable dir
  std::string dirname = get_environment_variable("NDUST_DATA");
  if(dirname.empty()) {
    std::cerr <<  "Environment variable NDUST_DATA not set" << '\n';
  }
  else {
    std::cerr <<  "Environment variable NDUST_DATA = " << dirname << '\n';
  }
  
  // Checks if option -o is called
  input_opt = get_cmd_option(argv, argv + argc, "-d");
  if(input_opt) {
    // sets new prefix filename
    dirname = std::string(input_opt);
    std::cerr << "\nUsing DATA dir from command line -d = " << dirname << '\n';
    if(dirname.empty()) {
      std::string dirname = "";
    }
  }
  
  std::cerr << "\n[ii] Grid file: " << grid_filename << '\n';
  std::cerr << "\n[ii] Plasma file: " << plasma_filename << '\n';
  std::cerr << "\n[ii] Nanoparticle file: " << nano_filename << '\n';
 
  // init logger
  blog::init(dirname+nano_filename);
   
  logging::add_common_attributes();

  src::severity_logger< severity_level > lg;
  BOOST_LOG_SEV(lg, info) << "Logging started for NEvo";
  
  BOOST_LOG_SEV(lg, info) << "Grid file: " << grid_filename;
  
  BOOST_LOG_SEV(lg, info) << "Output directory: " << dirname;
    // init logger
//   blog::close();
  // add /
//   dirname += "/";
  
  // get num threads
  int num_threads = omp_get_num_threads();
  BOOST_LOG_SEV(lg, info) << "Number of threads = " << num_threads;

  // get max threads
  int max_num_threads = omp_get_max_threads();
  BOOST_LOG_SEV(lg, info) <<  "Max number of threads = " << max_num_threads;
  
  input_opt = get_cmd_option(argv, argv + argc, "-t");  
  if(input_opt) {
    // get num threads from input
    num_threads = std::stoi(input_opt);
    // checks if num_threads > max_num_threads
    if(num_threads>max_num_threads) {
      num_threads = 1;
    }
    // set num_threads
    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
      #pragma omp master
      BOOST_LOG_SEV(lg, info) << "Set number of threads = " << num_threads;
    }
  }

  BOOST_LOG_SEV(lg, info) <<  "Get omp nested = " << omp_get_nested();
  BOOST_LOG_SEV(lg, info) <<  "Get omp cancellation = " << omp_get_cancellation();

  clock_t begin_nevo = std::clock();

  NEvo nevo(dirname, grid_filename, plasma_filename, nano_filename, lg);
  nevo.open();
  nevo.read();
  nevo.evolve();
  nevo.write();
  nevo.close();

  clock_t end_nevo = std::clock();
  double elapsed_secs = double(end_nevo - begin_nevo) / CLOCKS_PER_SEC;
  std::cout << "\n[ii] NEvo elapsed time: "
            << elapsed_secs << "\n\n";
  BOOST_LOG_SEV(lg, info)<< "NEvo elapsed time: " << elapsed_secs;

  return 0;
}
