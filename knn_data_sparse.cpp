//
// 
//                This source code is part of
// 
//                        M D S C T K
// 
//       Molecular Dynamics Spectral Clustering ToolKit
// 
//                        VERSION 1.1.2
// Written by Joshua L. Phillips.
// Copyright (c) 2013, Joshua L. Phillips.
// check out http://github.com/douradopalmares/mdsctk/ for more information.
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
// 
// If you want to redistribute modifications, please consider that
// derived work must not be called official MDSCTK. Details are found
// in the README & LICENSE files - if they are missing, get the
// official version at github.com/douradopalmares/mdsctk/.
// 
// To help us fund MDSCTK development, we humbly ask that you cite
// the papers on the package - you can find them in the top README file.
// 
// For more info, check our website at http://github.com/douradopalmares/mdsctk/
// 
//

// Standard
// C
#include <strings.h>
#include <stdlib.h>
#include <math.h>
// C++
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <ctime>

// Boost
#include <boost/program_options.hpp>

// OpenMP
#include <omp.h>

// Local
#include "config.h"
#include "mdsctk.h"

namespace po = boost::program_options;
using namespace std;

// You can change this function to utilize different
// distance metrics. The current implementation
// uses Euclidean distance...
double distance(int ref_size, int* ref_index, double* ref_data,
		int fit_size, int* fit_index, double* fit_data) {
  double value = 0.0;
  int ref = 0;
  int fit = 0;

  for (ref = 0; ref < ref_size; ref++) {
    while (fit < fit_size && fit_index[fit] < ref_index[ref]) {
      value += (fit_data[fit] * fit_data[fit]);
      fit++;
    }
    if (fit < fit_size && ref_index[ref] == fit_index[fit]) {
      value += ((ref_data[ref] - fit_data[fit]) *
		(ref_data[ref] - fit_data[fit]));
      fit++;
    }
    else {
      value += (ref_data[ref] * ref_data[ref]);
    }
  }
  for (;fit < fit_size;fit++) {
    value += (fit_data[fit] * fit_data[fit]);
  }
  return (sqrt(value));
}

vector<double> fits;

bool compare(int left, int right) {
  return fits[left] < fits[right];
}

int main(int argc, char* argv[]) {

  const char* program_name = "knn_data_sparse";
  bool optsOK = true;
  copyright(program_name);
  cout << "   Computes the k nearest neighbors of all pairs of" << endl;
  cout << "   vectors in the given sparse binary data files." << endl;
  cout << endl;
  cout << "   Use -h or --help to see the complete list of options." << endl;
  cout << endl;

  // Option vars...
  int nthreads = 0;
  int k = 0;
  string ref_index_filename;
  string ref_data_filename;
  string fit_index_filename;
  string fit_data_filename;
  string d_filename;
  string i_filename;

  // Declare the supported options.
  po::options_description cmdline_options;
  po::options_description program_options("Program options");
  program_options.add_options()
    ("help,h", "show this help message and exit")
    ("threads,t", po::value<int>(&nthreads)->default_value(2), "Input: Number of threads to start (int)")
    ("knn,k", po::value<int>(&k), "Input:  K-nearest neighbors (int)")
    ("reference-index-file,R", po::value<string>(&ref_index_filename)->default_value("reference.svi"), "Input:  Reference index file (string:filename)")
    ("reference-data-file,r", po::value<string>(&ref_data_filename)->default_value("reference.svd"), "Input:  Reference data file (string:filename)")
    ("fit-index-file,F", po::value<string>(&fit_index_filename), "Input:  Fitting index file (string:filename)")
    ("fit-data-file,f", po::value<string>(&fit_data_filename), "Input:  Fitting data file (string:filename)")
    ("distance-file,d", po::value<string>(&d_filename)->default_value("distances.dat"), "Output: K-nn distances file (string:filename)")
    ("index-file,i", po::value<string>(&i_filename)->default_value("indices.dat"), "Output: K-nn indices file (string:filename)")    
    ;
  cmdline_options.add(program_options);

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
  po::notify(vm);    

  if (vm.count("help")) {
    cout << "usage: " << program_name << " [options]" << endl;
    cout << cmdline_options << endl;
    return 1;
  }
  if (!vm.count("knn")) {
    cout << "ERROR: --knn not supplied." << endl;
    cout << endl;
    optsOK = false;
  }
  if (!vm.count("fit-index-file"))
    fit_index_filename = ref_index_filename;
  if (!vm.count("fit-data-file"))
    fit_data_filename = ref_data_filename;

  if (!optsOK) {
    return -1;
  }

  cout << "Running with the following options:" << endl;
  cout << "threads =              " << nthreads << endl;
  cout << "knn =                  " << k << endl;
  cout << "reference-index-file = " << ref_index_filename << endl;
  cout << "reference-file =       " << ref_data_filename << endl;
  cout << "fit-index-file =       " << fit_index_filename << endl;
  cout << "fit-file =             " << fit_data_filename << endl;
  cout << "distance-file =        " << d_filename << endl;
  cout << "index-file =           " << i_filename << endl;
  cout << endl;
  
  // Local vars
  vector<int> ref_size;
  vector<int> fit_size;
  vector<int*> ref_index;
  vector<int*> fit_index;
  vector<double*> ref_data;
  vector<double*> fit_data;
  int k1 = k + 1;
  int update_interval = 1;
  vector<double> keepers;
  vector<int> permutation;
  ifstream index;
  ifstream data;
  ofstream distances;
  ofstream indices;
  int n;
  int *myindex;
  double *mydata;

  // Setup threads
  omp_set_num_threads(nthreads);

  // Read coordinates
  cout << "Reading reference coordinates from files...";
  index.open(ref_index_filename.c_str());
  data.open(ref_data_filename.c_str());
  index.read((char*) &n, sizeof(int) / sizeof(char));
  while (!index.eof()) {
    myindex = new int[n];
    mydata = new double[n];
    index.read((char*) myindex, (sizeof(int) / sizeof(char)) * n);
    data.read((char*) mydata, (sizeof(double) / sizeof(char)) * n);
    ref_size.push_back(n);
    ref_index.push_back(myindex);
    ref_data.push_back(mydata);
    index.read((char*) &n, sizeof(int) / sizeof(char));
  }
  index.close();
  data.close();  
  cout << "done." << endl;

  cout << "Reading fitting coordinates from files...";
  index.open(fit_index_filename.c_str());
  data.open(fit_data_filename.c_str());
  index.read((char*) &n, sizeof(int) / sizeof(char));
  while (!index.eof()) {
    myindex = new int[n];
    mydata = new double[n];
    index.read((char*) myindex, (sizeof(int) / sizeof(char)) * n);
    data.read((char*) mydata, (sizeof(double) / sizeof(char)) * n);
    fit_size.push_back(n);
    fit_index.push_back(myindex);
    fit_data.push_back(mydata);
    index.read((char*) &n, sizeof(int) / sizeof(char));
  }
  index.close();
  data.close();  
  cout << "done." << endl;

  // Open output files
  distances.open(d_filename.c_str());
  indices.open(i_filename.c_str());

  // Allocate vectors for storing the RMSDs for a structure
  fits.resize(ref_index.size());
  permutation.resize(ref_index.size());

  // Fix k if number of frames is too small
  if (ref_index.size()-1 < k)
    k = ref_index.size()-1;
  k1 = k + 1;
  keepers.resize(k1);

  // Timers...
  time_t start = std::time(0);
  time_t last = start;

  // Compute fits
  for (int fit_frame = 0; fit_frame < fit_index.size(); fit_frame++) {
    
    // Update user of progress
    if (std::time(0) - last > update_interval) {
      last = std::time(0);
      time_t eta = start + ((last-start) * fit_index.size() / fit_frame);
      cout << "\rFrame: " << fit_frame << ", will finish " 
	   << string(std::ctime(&eta)).substr(0,20);
      cout.flush();
    }

    // Do Work
#pragma omp parallel for
    for (int ref_frame = 0; ref_frame < ref_index.size(); ref_frame++)
      fits[ref_frame] =
	distance(ref_size[ref_frame],ref_index[ref_frame],ref_data[ref_frame],
		 fit_size[fit_frame],fit_index[fit_frame],fit_data[fit_frame]);

    // Sort
    int x = 0;
    for (vector<int>::iterator p_itr = permutation.begin();
	 p_itr != permutation.end(); p_itr++)
      (*p_itr) = x++;    
    partial_sort(permutation.begin(), permutation.begin()+k1,
		 permutation.end(), compare);
    for (int x = 0; x < k1; x++)
      keepers[x] = (double) fits[permutation[x]];

    // Write out closest k RMSD alignment scores and indices
    distances.write((char*) &(keepers[1]), (sizeof(double) / sizeof(char)) * k);
    indices.write((char*) &(permutation[1]), (sizeof(int) / sizeof(char)) * k);
  }

  cout << endl << endl;

  // Clean coordinates
  for (vector<int*>::iterator itr = ref_index.begin();
       itr != ref_index.end(); itr++) delete [] (*itr);
  for (vector<int*>::iterator itr = fit_index.begin();
       itr != fit_index.end(); itr++) delete [] (*itr);
  for (vector<double*>::iterator itr = ref_data.begin();
       itr != ref_data.end(); itr++) delete [] (*itr);
  for (vector<double*>::iterator itr = fit_data.begin();
       itr != fit_data.end(); itr++) delete [] (*itr);

  return 0;
}
