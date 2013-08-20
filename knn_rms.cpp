//
// 
//                This source code is part of
// 
//                        M D S C T K
// 
//       Molecular Dynamics Spectral Clustering ToolKit
// 
//                        VERSION 1.1.1
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

// GROMACS
#include <gromacs/tpxio.h>
#include <gromacs/xtcio.h>
#include <gromacs/do_fit.h>

// Local
#include "config.h"

using namespace std;
typedef real (*coord_array)[3];

int main(int argc, char* argv[]) {
  
  if (argc != 6) {
    cerr << endl;
    cerr << "   MDSCTK " << MDSCTK_VERSION_MAJOR << "." << MDSCTK_VERSION_MINOR << endl;
    cerr << "   Copyright (C) 2013 Joshua L. Phillips" << endl;
    cerr << "   MDSCTK comes with ABSOLUTELY NO WARRANTY; see LICENSE for details." << endl;
    cerr << "   This is free software, and you are welcome to redistribute it" << endl;
    cerr << "   under certain conditions; see README.md for details." << endl;
    cerr << endl;
    cerr << "Usage: " << argv[0] << " [# threads] [k] [topology file] [reference xtc file] [fitting xtc file]" << endl;
    cerr << "   Computes the k nearest neighbors of all reference" << endl;
    cerr << "   structures in the given xtc file for each structure" << endl;
    cerr << "   in the given fitting xtc file. (Use the same file for" << endl;
    cerr << "   fitting and reference to make a symmetric comparison.)" << endl;
    cerr << "   A topology PDB file should be provided for determining" << endl;
    cerr << "   the mass of each atom." << endl;
    cerr << endl;
    return -1;
  }

  // Local vars
  int step = 1;
  float time = 0.0;
  matrix box;
  float prec = 0.001;
  char buf[256];
  t_topology top;
  int ePBC;
  int nthreads = atoi(argv[1]);
  int natoms = 0;
  int k = atoi(argv[2]);
  int k1 = k + 1;
  int update_interval = 1;
  const char* top_filename = argv[3];
  const char* ref_filename = argv[4];
  const char* fit_filename = argv[5];
  t_fileio *ref_file;
  t_fileio *fit_file;
  rvec *mycoords = NULL;
  gmx_bool bOK = 1;
  gsl_vector *keepers = NULL;
  gsl_permutation *permutation = NULL;
  ofstream distances;
  ofstream indices;
  vector<coord_array> *ref_coords = NULL;
  vector<coord_array> *fit_coords = NULL;
  real *weights = NULL;
  gsl_vector *fits = NULL;
  int *fit_indices = NULL;

  // Setup threads
  omp_set_num_threads(nthreads);

  // Get number of atoms and initialize weights
  cerr << "Reading topology information from " << top_filename << " ... ";
  read_tps_conf(top_filename, buf, &top, &ePBC, &mycoords,
		NULL, box, TRUE);
  cerr << "done." << endl;
  delete [] mycoords;

  ref_file = open_xtc(ref_filename,"r");
  read_first_xtc(ref_file,&natoms, &step, &time, box, &mycoords, &prec, &bOK);
  close_xtc(ref_file);
  if (natoms != top.atoms.nr) {
    cerr << "*** ERROR ***" << endl;
    cerr << "Number of atoms in topology file ("
	 << top.atoms.nr << ") "
	 << "does not match the number of atoms "
	 << "in the XTC file (" << ref_filename << " : " << natoms << ")."
	 << endl;
    exit(4);
  }
  fit_file = open_xtc(fit_filename,"r");
  read_first_xtc(fit_file,&natoms, &step, &time, box, &mycoords, &prec, &bOK);
  close_xtc(fit_file);
  if (natoms != top.atoms.nr) {
    cerr << "*** ERROR ***" << endl;
    cerr << "Number of atoms in topology file ("
	 << top.atoms.nr << ") "
	 << "does not match the number of atoms "
	 << "in the XTC file (" << fit_filename << " : " << natoms << ")."
	 << endl;
    exit(4);
  }
  ref_coords = new vector<coord_array>;
  fit_coords = new vector<coord_array>;
  weights = new real[natoms];
  for (int x = 0; x < natoms; x++) weights[x] = top.atoms.atom[x].m;

  // Read coordinates and weight-center all structures
  cerr << "Reading reference coordinates from file: " << ref_filename << " ... ";
  ref_file = open_xtc(ref_filename,"r");
  mycoords = new rvec[natoms];
  while (read_next_xtc(ref_file, natoms, &step, &time, box, mycoords, &prec, &bOK)) {
    reset_x(natoms,NULL,natoms,NULL,mycoords,weights);
    ref_coords->push_back(mycoords);
    mycoords = new rvec[natoms];
  }
  close_xtc(ref_file);
  cerr << "done." << endl;

  cerr << "Reading fitting coordinates from file: " << fit_filename << " ... ";
  fit_file = open_xtc(fit_filename,"r");
  while (read_next_xtc(fit_file, natoms, &step, &time, box, mycoords, &prec, &bOK)) {
    reset_x(natoms,NULL,natoms,NULL,mycoords,weights);
    fit_coords->push_back(mycoords);
    mycoords = new rvec[natoms];
  }
  close_xtc(fit_file);
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
    for (int ref_frame = 0; ref_frame < ref_coords->size(); ref_frame++) {
      do_fit(natoms,weights,
	     (*fit_coords)[fit_frame],
	     (*ref_coords)[ref_frame]);
      gsl_vector_set(fits,ref_frame,rmsdev(natoms,weights,
					   (*fit_coords)[fit_frame],
					   (*ref_coords)[ref_frame]) * 10.0);
    }

    // Sort
    gsl_permutation_init(permutation);
    gsl_sort_vector_index(permutation,fits);
    for (int x = 0; x < k1; x++) {
      gsl_vector_set(keepers,x,gsl_vector_get(fits,gsl_permutation_get(permutation,x)));
      fit_indices[x] = (int) gsl_permutation_get(permutation,x);
    }
    // Write out closest k RMSD alignment scores and indices
    distances.write((char*) &(keepers->data[1]), (sizeof(double)/sizeof(char)) * k);
    indices.write((char*) &(fit_indices[1]), (sizeof(int)/sizeof(char)) * k);

  }

  cerr << "\rWorking: " << 100.0 << "%" << endl;

  // Clean coordinates
  for (vector<coord_array>::iterator itr = ref_coords->begin();
       itr != ref_coords->end(); itr++) delete [] (*itr);
  for (vector<coord_array>::iterator itr = fit_coords->begin();
       itr != fit_coords->end(); itr++) delete [] (*itr);
  delete ref_coords;
  delete fit_coords;
  delete [] weights;
  delete [] fit_indices;

  gsl_vector_free(fits);
  gsl_vector_free(keepers);
  gsl_permutation_free(permutation);

  return 0;
}
