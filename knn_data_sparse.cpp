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

// GSL Tools
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort_vector.h>

// Standard
#include <iostream>
#include <fstream>
#include <vector>
#include <strings.h>
#include <stdlib.h>
#include <math.h>

// OpenMP
#include <omp.h>

// Local
#include "config.h"

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
      // cout << "F1: " 
      // 	   << ref_index[ref] << " " << fit_index[fit] << " "
      // 	   << ref_data[ref] << " " << fit_data[fit] << endl;
      fit++;
    }
    if (fit < fit_size && ref_index[ref] == fit_index[fit]) {
      value += ((ref_data[ref] - fit_data[fit]) *
		(ref_data[ref] - fit_data[fit]));
      // cout << "M:  " 
      // 	   << ref_index[ref] << " " << fit_index[fit] << " "
      // 	   << ref_data[ref] << " " << fit_data[fit] << endl;
      fit++;
    }
    else {
      // cout << "R:  " 
      // 	   << ref_index[ref] << " " << fit_index[fit] << " "
      // 	   << ref_data[ref] << " " << fit_data[fit] << endl;
      value += (ref_data[ref] * ref_data[ref]);
    }
  }
  for (;fit < fit_size;fit++) {
      // cout << "F2: " 
      // 	   << ref_index[ref] << " " << fit_index[fit] << " "
      // 	   << ref_data[ref] << " " << fit_data[fit] << endl;
    value += (fit_data[fit] * fit_data[fit]);
  }
  return (sqrt(value));
}

int main(int argc, char* argv[]) {
  
  if (argc != 3) {
    cerr << endl;
    cerr << "   MDSCTK " << MDSCTK_VERSION_MAJOR << "." << MDSCTK_VERSION_MINOR << endl;
    cerr << endl;
    cerr << "Usage: " << argv[0] << " [# threads] [k]" << endl;
    cerr << "   Computes the k nearest neighbors of all pairs of" << endl;
    cerr << "   vectors in the given sparse binary data files." << endl;
    cerr << endl;
    return -1;
  }

  // Local vars
  int nthreads = 0;
  vector<int> ref_size;
  vector<int> fit_size;
  vector<int*> ref_index;
  vector<int*> fit_index;
  vector<double*> ref_data;
  vector<double*> fit_data;
  gsl_vector *fits = NULL;
  int *fit_indices = NULL;
  nthreads = atoi(argv[1])-1;
  int k = atoi(argv[2]);
  int k1 = k + 1;
  int update_interval = 1;
  gsl_vector *keepers = NULL;
  gsl_permutation *permutation = NULL;
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
  cerr << "Reading reference coordinates from files...";
  index.open("index.dat");
  data.open("data.dat");
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
  cerr << "done." << endl;

  cerr << "Reading fitting coordinates from files...";
  index.open("index.dat");
  data.open("data.dat");
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
  cerr << "done." << endl;

  // Open output files
  distances.open("distances.dat");
  indices.open("indices.dat");

  // Allocate vectors for storing the RMSDs for a structure
  fits = gsl_vector_calloc(ref_index.size());
  permutation = gsl_permutation_calloc(ref_index.size());
  fit_indices = new int[ref_index.size()];

  // Fix k if number of frames is too small
  if (ref_index.size()-1 < k)
    k = ref_index.size()-1;
  k1 = k + 1;
  keepers = gsl_vector_calloc(k1);

  // Get update frequency
  cerr.precision(8);
  cerr.setf(ios::fixed,ios::floatfield);
  update_interval = ceil(sqrt(fit_index.size()));

  // // Test
  // cout << endl;
  // cout << "Testing..." << endl;
  // cout << "Reference: " << ref_size[0] << endl;
  // cout << "Ref Ind:" << endl;
  // for (int i = 0; i < ref_size[0]; i++)
  //   cout << ref_index[0][i] << " ";
  // cout << endl;
  // cout << "Ref Data:" << endl;
  // for (int i = 0; i < ref_size[0]; i++)
  //   cout << ref_data[0][i] << " ";
  // cout << endl;
  // cout << "Fitting: " << fit_size[1] << endl;
  // cout << "Fit Ind:" << endl;
  // for (int i = 0; i < fit_size[1]; i++)
  //   cout << fit_index[1][i] << " ";
  // cout << endl;
  // cout << "Fit Data:" << endl;
  // for (int i = 0; i < fit_size[1]; i++)
  //   cout << fit_data[1][i] << " ";
  // cout << endl;
  // cout << "Distance:" << distance(ref_size[1],ref_index[1],ref_data[1],
  // 				  fit_size[0],fit_index[0],fit_data[0]) << endl;
  // cout << endl;

  // return 0;

  // Compute fits
  for (int fit_frame = 0; fit_frame < fit_index.size(); fit_frame++) {
    
    // Update user of progress
    if (fit_frame % update_interval == 0) {
      cerr << "\rWorking: " << (((double) fit_frame) / ((double) fit_index.size())) * 100.0 << "%";
      cerr.flush();
    }

    // Do Work
#pragma omp parallel for
    for (int ref_frame = 0; ref_frame < ref_index.size(); ref_frame++)
      gsl_vector_set(fits,ref_frame,
		     distance(ref_size[ref_frame],ref_index[ref_frame],ref_data[ref_frame],
			      fit_size[fit_frame],fit_index[fit_frame],fit_data[fit_frame]));

    // Sort
    gsl_permutation_init(permutation);
    gsl_sort_vector_index(permutation,fits);
    for (int x = 0; x < k1; x++) {
      gsl_vector_set(keepers,x,gsl_vector_get(fits,gsl_permutation_get(permutation,x)));
      fit_indices[x] = (int) gsl_permutation_get(permutation,x);
    }
    // Write out closest k RMSD alignment scores and indices
    distances.write((char*) &(keepers->data[1]), (sizeof(double) / sizeof(char)) * k);
    indices.write((char*) &(fit_indices[1]), (sizeof(int) / sizeof(char)) * k);
  }

  cerr << "\rWorking: " << 100.0 << "%" << endl;

  // Clean coordinates
  for (vector<int*>::iterator itr = ref_index.begin();
       itr != ref_index.end(); itr++) delete [] (*itr);
  for (vector<int*>::iterator itr = fit_index.begin();
       itr != fit_index.end(); itr++) delete [] (*itr);
  for (vector<double*>::iterator itr = ref_data.begin();
       itr != ref_data.end(); itr++) delete [] (*itr);
  for (vector<double*>::iterator itr = fit_data.begin();
       itr != fit_data.end(); itr++) delete [] (*itr);
  delete [] fit_indices;

  gsl_vector_free(fits);
  gsl_vector_free(keepers);
  gsl_permutation_free(permutation);

  return 0;
}
