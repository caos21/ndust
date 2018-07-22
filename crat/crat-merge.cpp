#include <iostream>
#include <string>
#include <fstream>

// OpenMP
#include <omp.h>

// io functions
#include "io.h"

// logger
#include "../include/log.h"

// class CRate
#include "../include/CRate.h"

int main(int argc, char **argv) {
  // Prints header
  print_header(std::cout);

  // Print help
  show_help_merge(argc, argv);

  // Checks if option -o is called
  char *input_opt;
  input_opt = get_cmd_option(argv, argv + argc, "-o");
  std::string prefix_filename = "";
  if(input_opt) {
    // sets new prefix filename
    prefix_filename.append(input_opt);
  }
  else {
    std::cerr << "\n[ee] No prefix given. Terminate." << '\n';
    return -1;
  }
  
  // Checks if option -i is called
  input_opt = get_cmd_option(argv, argv + argc, "-i");
  std::string input_filename = "";
  if(input_opt) {
    // sets new prefix filename
    input_filename.append(input_opt);
  }
  else {
    std::cerr << "\n[ww] No input list given.\n";
    input_filename = std::string(prefix_filename + ".list");
    std::cerr << "\n[ww] Using input list : " << input_filename << '\n';
  }


  // Check environment variable dir
  std::string dirname = get_environment_variable("NDUST_DATA");
  if(dirname.empty()) {
    std::cerr << "\n[ii] Environment variable NDUST_DATA not set" << '\n';
  }
  else {
    std::cerr << "\n[ii] Environment variable NDUST_DATA = " << dirname << '\n';
  }

  // Checks if option -d is called
  input_opt = get_cmd_option(argv, argv + argc, "-d");
  if(input_opt) {
    // sets new prefix filename
    dirname = std::string(input_opt);
    std::cerr << "\n[ii] Using DATA dir from command line -d = " << dirname << '\n';
    if(dirname.empty()) {
      std::string dirname = "";
    }
  }

  std::cerr << "\n[ii] Input file list: " << input_filename << '\n';

  std::cerr << "\n[ii] Prefix for output files: " << prefix_filename << '\n';
  
  // init logger
  blog::init(std::string(dirname+prefix_filename+"-merge"));

  logging::add_common_attributes();

  src::severity_logger< severity_level > lg;
  BOOST_LOG_SEV(lg, info) << "Logging started for CRat merge";

  BOOST_LOG_SEV(lg, info) << "Input file list: " << input_filename;

  BOOST_LOG_SEV(lg, info) << "Prefix for output files: " << prefix_filename;
    
  BOOST_LOG_SEV(lg, info) << "Output directory: " << dirname;

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
      #pragma omp single
      BOOST_LOG_SEV(lg, info) << "Set number of threads = " << num_threads;
    }
  }

  BOOST_LOG_SEV(lg, info) <<  "Get omp nested = " << omp_get_nested();
  BOOST_LOG_SEV(lg, info) <<  "Get omp cancellation = " << omp_get_cancellation();

  clock_t begin_crate = std::clock();


  // get list of h5 files
  std::fstream input_filelist(dirname+input_filename, std::fstream::in);
  std::vector<std::string> sfilelist;
  if(input_filelist.is_open()) {
    BOOST_LOG_SEV(lg, info) <<  "Processing input file list: " << input_filename;
    std::string sfilename;
    while (input_filelist >> sfilename) {
      BOOST_LOG_SEV(lg, info) <<  "File : " << sfilename;
      sfilelist.push_back(sfilename);
    }
    if (sfilelist.size()==0) {
      std::cerr << "\n[ee] No elements in input file list: " << input_filename << '\n';
      BOOST_LOG_SEV(lg, fatal) <<  "No elements in input file list error: " << input_filename;
      std::cerr << "\n[ee] Terminate." << '\n';
      return -1;
    }
  }
  else {
    std::cerr << "\n[ee] Input file list error: " << input_filename << '\n';
    BOOST_LOG_SEV(lg, fatal) <<  "Input file list error: " << input_filename;
    std::cerr << "\n[ee] Terminate." << '\n';
    return -1;
  }
  
  CRate crate(dirname, prefix_filename, lg);
  crate.open();
  crate.read();
  crate.compute_list(sfilelist);
  crate.write();
  crate.close();

  clock_t end_crate = std::clock();
  double elapsed_secs = double(end_crate - begin_crate) / CLOCKS_PER_SEC;
  std::cout << "\n[ii] CRat elapsed time: "
            << elapsed_secs << '\n';
  BOOST_LOG_SEV(lg, info)<< "CRat reader elapsed time: " << elapsed_secs;

  return 0;
}
