//
// 
//                This source code is part of
// 
//                        M D S C T K
// 
//       Molecular Dynamics Spectral Clustering ToolKit
// 
//                        VERSION 1.2.2
// Written by Joshua L. Phillips.
// Copyright (c) 2012-2014, Joshua L. Phillips.
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

// Local
#include "config.h"
#include "mdsctk.h"

// You can change this function to utilize different
// distance metrics. The current implementation
// uses Euclidean distance...
double (*distance)(int ref_size, int* ref_index, double* ref_data,
		   int fit_size, int* fit_index, double* fit_data) =
  euclidean_distance_sparse;

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
  int blksize = 128;
  string ref_index_filename;
  string ref_data_filename;
  string fit_index_filename;
  string fit_data_filename;
  string d_filename;
  string i_filename;
  bool sort;

  // Declare the supported options.
  po::options_description cmdline_options;
  po::options_description program_options("Program options");
  program_options.add_options()
    ("help,h", "show this help message and exit")
    ("threads,t", po::value<int>(&nthreads)->default_value(omp_get_max_threads()>omp_get_num_procs()?omp_get_num_procs():omp_get_max_threads()), "Input:  Number of threads to start (int)")
    ("knn,k", po::value<int>(&k), "Input:  K-nearest neighbors (int)")
    ("sort,s",po::value<bool>(&sort)->default_value(true),"Input:  Find K-nn,false=full distance matix (bool)")
    ("block-size,b", po::value<int>(&blksize)->default_value(128), "Input:  Workgroup block size in # frames (int)")
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
  if (!vm.count("knn") && sort) {
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
  cout << "sort =                 " << sort << endl;
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
  permutation<double> *fits;
  ifstream index;
  ifstream data;
  ofstream distances;
  ofstream indices;
  int n;
  int *myindex;
  double *mydata;
  int max_blks = 1024 * 1024 * 1024 / 8;

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
  max_blks /= fit_index.size();
  if (blksize > max_blks)
    blksize = max_blks;
  if (blksize > fit_index.size())
    blksize = fit_index.size();
  cout << "Block size: " << blksize << endl;     
  fits = new permutation<double>[blksize];
  for (int x = 0; x < blksize; x++)
    fits[x].data.resize(ref_index.size());

  // Fix k if number of frames is too small
  if (ref_index.size()-1 < k)
    k = ref_index.size()-1;
  k1 = k + 1;

  // Timers...
  time_t start = std::time(0);
  time_t last = start;

  // Compute fits
  int remainder = fit_index.size() % blksize;
#pragma omp parallel for
  for (int frame = 0; frame < remainder; frame++) {
    for (int ref_frame = 0; ref_frame < ref_index.size(); ref_frame++)
      fits[frame].data[ref_frame] =
	::distance(ref_size[ref_frame],ref_index[ref_frame],ref_data[ref_frame],
		   fit_size[frame],fit_index[frame],fit_data[frame]);
    if (sort)
      fits[frame].sort(k1);
  }

  // Write out closest k RMSD alignment scores and indices
  for (int frame = 0; frame < remainder; frame++) {
    if (sort) {
      distances.write((char*) &(fits[frame].data[1]), (sizeof(double)/sizeof(char)) * k);
      indices.write((char*) &(fits[frame].indices[1]), (sizeof(int)/sizeof(char)) * k);
    }
    else {
      distances.write((char*) &(fits[frame].data[0]), (sizeof(double)/sizeof(char)) * ref_index.size());
      indices.write((char*) &(fits[frame].indices[0]), (sizeof(int)/sizeof(char)) * ref_index.size());
    }
  }
  for (int fit_frame = remainder; fit_frame < fit_index.size(); fit_frame += blksize) {
    
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
    for (int frame = 0; frame < blksize; frame++) {
      for (int ref_frame = 0; ref_frame < ref_index.size(); ref_frame++)
	fits[frame].data[ref_frame] =
	  ::distance(ref_size[ref_frame],ref_index[ref_frame],ref_data[ref_frame],
		     fit_size[frame+fit_frame],fit_index[frame+fit_frame],fit_data[frame+fit_frame]);
      if (sort)
	fits[frame].sort(k1);
    }

    // Write out closest k RMSD alignment scores and indices
    for (int frame = 0; frame < blksize; frame++) {
      if (sort) {
	distances.write((char*) &(fits[frame].data[1]), (sizeof(double)/sizeof(char)) * k);
	indices.write((char*) &(fits[frame].indices[1]), (sizeof(int)/sizeof(char)) * k);
      }
      else {
	distances.write((char*) &(fits[frame].data[0]), (sizeof(double)/sizeof(char)) * ref_index.size());
	indices.write((char*) &(fits[frame].indices[0]), (sizeof(int)/sizeof(char)) * ref_index.size());
      }
    }
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
