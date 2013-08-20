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
double distance(int size, double* reference, double* fitting) {
  double value = 0.0;
  for (int x = 0; x < size; x++)
    value += (reference[x] - fitting[x]) * (reference[x] - fitting[x]);
  return (sqrt(value));
}

int main(int argc, char* argv[]) {
  
  if (argc != 6) {
    cerr << endl;
    cerr << "   MDSCTK " << MDSCTK_VERSION_MAJOR << "." << MDSCTK_VERSION_MINOR << endl;
    cerr << "   Copyright (C) 2013 Joshua L. Phillips" << endl;
    cerr << "   MDSCTK comes with ABSOLUTELY NO WARRANTY; see LICENSE for details." << endl;
    cerr << "   This is free software, and you are welcome to redistribute it" << endl;
    cerr << "   under certain conditions; see README.md for details." << endl;
    cerr << endl;
    cerr << "Usage: " << argv[0] << " [# threads] [k] [vector size] [reference data file] [fitting data file]" << endl;
    cerr << "   Computes the k nearest neighbors of all pairs of" << endl;
    cerr << "   vectors in the given binary data files." << endl;
    cerr << endl;
    return -1;
  }

  // Local vars
  int nthreads = 0;
  int vector_size = 0;
  vector<double*> *ref_coords = NULL;
  vector<double*> *fit_coords = NULL;
  gsl_vector *fits = NULL;
  int *fit_indices = NULL;
  nthreads = atoi(argv[1])-1;
  int k = atoi(argv[2]);
  int k1 = k + 1;
  vector_size = atoi(argv[3]);
  int update_interval = 1;
  const char* ref_file = argv[4];
  const char* fit_file = argv[5];
  gsl_vector *keepers = NULL;
  gsl_permutation *permutation = NULL;
  ofstream distances;
  ofstream indices;

  // Setup threads
  omp_set_num_threads(nthreads);

  ref_coords = new vector<double*>;
  fit_coords = new vector<double*>;

  // Read coordinates
  cerr << "Reading reference coordinates from file: " << ref_file << " ... ";
  ifstream myfile;
  myfile.open(ref_file);
  double* mycoords = new double[vector_size];
  myfile.read((char*) mycoords, sizeof(double) * vector_size);
  while (!myfile.eof()) {
    ref_coords->push_back(mycoords);
    mycoords = new double[vector_size];
    myfile.read((char*) mycoords, sizeof(double) * vector_size);
  }
  myfile.close();
  cerr << "done." << endl;

  cerr << "Reading fitting coordinates from file: " << fit_file << " ... ";
  myfile.open(fit_file);
  myfile.read((char*) mycoords, sizeof(double) * vector_size);
  while (!myfile.eof()) {
    fit_coords->push_back(mycoords);
    mycoords = new double[vector_size];
    myfile.read((char*) mycoords, sizeof(double) * vector_size);
  }
  myfile.close();
  delete [] mycoords;
  mycoords = NULL;
  cerr << "done." << endl;

  // Open output files
  distances.open("distances.dat");
  indices.open("indices.dat");

  // Allocate vectors for storing the RMSDs for a structure
  fits = gsl_vector_calloc(ref_coords->size());
  permutation = gsl_permutation_calloc(ref_coords->size());
  fit_indices = new int[ref_coords->size()];

  // Fix k if number of frames is too small
  if (ref_coords->size()-1 < k)
    k = ref_coords->size()-1;
  k1 = k + 1;
  keepers = gsl_vector_calloc(k1);

  // Get update frequency
  // int update_interval = (int) floor(sqrt((double) coords.size()));
  cerr.precision(8);
  cerr.setf(ios::fixed,ios::floatfield);
  update_interval = ceil(sqrt(fit_coords->size()));

  // Compute fits
  for (int fit_frame = 0; fit_frame < fit_coords->size(); fit_frame++) {
    
    // Update user of progress
    if (fit_frame % update_interval == 0) {
      cerr << "\rWorking: " << (((double) fit_frame) / ((double) fit_coords->size())) * 100.0 << "%";
      cerr.flush();
    }

    // Do Work
#pragma omp parallel for
    for (int ref_frame = 0; ref_frame < ref_coords->size(); ref_frame++)
      gsl_vector_set(fits,ref_frame,distance(vector_size,(*fit_coords)[fit_frame],
					     (*ref_coords)[ref_frame]));


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
  for (vector<double*>::iterator itr = ref_coords->begin();
       itr != ref_coords->end(); itr++) delete [] (*itr);
  for (vector<double*>::iterator itr = fit_coords->begin();
       itr != fit_coords->end(); itr++) delete [] (*itr);
  delete ref_coords;
  delete fit_coords;
  delete [] fit_indices;

  gsl_vector_free(fits);
  gsl_vector_free(keepers);
  gsl_permutation_free(permutation);

  return 0;
}
