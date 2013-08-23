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
#include <ctime>

// Boost
#include <boost/program_options.hpp>

// OpenMP
#include <omp.h>

// GROMACS
#include <gromacs/tpxio.h>
#include <gromacs/xtcio.h>
#include <gromacs/do_fit.h>

// Local
#include "config.h"
#include "mdsctk.h"

typedef real (*coord_array)[3];

namespace po = boost::program_options;
using namespace std;

vector<real> fits;

bool compare(int left, int right) {
  return fits[left] < fits[right];
}

int main(int argc, char* argv[]) {

  const char* program_name = "knn_rms";
  bool optsOK = true;
  copyright(program_name);
  cout << "   Computes the k nearest neighbors of all reference" << endl;
  cout << "   structures in the given xtc file for each structure" << endl;
  cout << "   in the given fitting xtc file. (Uses the same file for" << endl;
  cout << "   as reference by default to make a symmetric comparison.)" << endl;
  cout << "   A topology PDB file should be provided for determining" << endl;
  cout << "   the mass of each atom." << endl;
  cout << endl;
  cout << "   Use -h or --help to see the complete list of options." << endl;
  cout << endl;

  // Option vars...
  int nthreads = 0;
  int k = 0;
  string top_filename;
  string ref_filename;
  string fit_filename;
  string d_filename;
  string i_filename;

  // Declare the supported options.
  po::options_description cmdline_options;
  po::options_description program_options("Program options");
  program_options.add_options()
    ("help,h", "show this help message and exit")
    ("threads,t", po::value<int>(&nthreads)->default_value(2), "Input: Number of threads to start (int)")
    ("knn,k", po::value<int>(&k), "Input:  K-nearest neighbors (int)")
    ("topology-file,p", po::value<string>(&top_filename)->default_value("topology.pdb"), "Input:  Topology file [.pdb,.gro,.tpr] (string:filename)")
    ("reference-file,r", po::value<string>(&ref_filename)->default_value("reference.xtc"), "Input:  Reference [.xtc] file (string:filename)")
    ("fit-file,f", po::value<string>(&fit_filename), "Input:  Fitting [.xtc] file (string:filename)")
    ("distance-file,d", po::value<string>(&d_filename)->default_value("distances.dat"), "Output: K-nn distances file (string:filename)")
    ("index-file,i", po::value<string>(&i_filename)->default_value("indices.dat"), "Output: K-nn indices file (string:filename)")    
    ;
  cmdline_options.add(program_options);

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
  po::notify(vm);    

  if (vm.count("help")) {
    cout << "usage: " << program_name << " [options]" << endl;
    cout << endl;
    cout << cmdline_options << endl;
    return 1;
  }
  if (!vm.count("knn")) {
    cout << "ERROR: --knn not supplied." << endl;
    cout << endl;
    optsOK = false;
  }
  if (!vm.count("fit-file"))
    fit_filename = ref_filename;

  if (!optsOK) {
    return -1;
  }

  cout << "Running with the following options:" << endl;
  cout << "threads =        " << nthreads << endl;
  cout << "knn =            " << k << endl;
  cout << "topology-file =  " << top_filename << endl;
  cout << "reference-file = " << ref_filename << endl;
  cout << "fit-file =       " << fit_filename << endl;
  cout << "distance-file =  " << d_filename << endl;
  cout << "index-file =     " << i_filename << endl;
  cout << endl;
  
  // Local vars
  int step = 1;
  float time = 0.0;
  matrix box;
  float prec = 0.001;
  char buf[256];
  t_topology top;
  int ePBC;
  int natoms = 0;
  int k1 = k + 1;
  int update_interval = 1;
  t_fileio *ref_file;
  t_fileio *fit_file;
  rvec *mycoords = NULL;
  gmx_bool bOK = 1;
  vector<double> keepers;
  vector<int> permutation;
  ofstream distances;
  ofstream indices;
  vector<coord_array> *ref_coords = NULL;
  vector<coord_array> *fit_coords = NULL;
  real *weights = NULL;

  // Setup threads
  omp_set_num_threads(nthreads);

  // Get number of atoms and initialize weights
  cout << "Reading topology information from " << top_filename << " ... ";
  read_tps_conf(top_filename.c_str(), buf, &top, &ePBC, &mycoords,
		NULL, box, TRUE);
  cout << "done." << endl;
  delete [] mycoords;

  ref_file = open_xtc(ref_filename.c_str(),"r");
  read_first_xtc(ref_file,&natoms, &step, &time, box, &mycoords, &prec, &bOK);
  close_xtc(ref_file);
  if (natoms != top.atoms.nr) {
    cout << "*** ERROR ***" << endl;
    cout << "Number of atoms in topology file ("
	 << top.atoms.nr << ") "
	 << "does not match the number of atoms "
	 << "in the XTC file (" << ref_filename << " : " << natoms << ")."
	 << endl;
    exit(4);
  }
  fit_file = open_xtc(fit_filename.c_str(),"r");
  read_first_xtc(fit_file,&natoms, &step, &time, box, &mycoords, &prec, &bOK);
  close_xtc(fit_file);
  if (natoms != top.atoms.nr) {
    cout << "*** ERROR ***" << endl;
    cout << "Number of atoms in topology file ("
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
  cout << "Reading reference coordinates from file: " << ref_filename << " ... ";
  ref_file = open_xtc(ref_filename.c_str(),"r");
  mycoords = new rvec[natoms];
  while (read_next_xtc(ref_file, natoms, &step, &time, box, mycoords, &prec, &bOK)) {
    reset_x(natoms,NULL,natoms,NULL,mycoords,weights);
    ref_coords->push_back(mycoords);
    mycoords = new rvec[natoms];
  }
  close_xtc(ref_file);
  cout << "done." << endl;

  cout << "Reading fitting coordinates from file: " << fit_filename << " ... ";
  fit_file = open_xtc(fit_filename.c_str(),"r");
  while (read_next_xtc(fit_file, natoms, &step, &time, box, mycoords, &prec, &bOK)) {
    reset_x(natoms,NULL,natoms,NULL,mycoords,weights);
    fit_coords->push_back(mycoords);
    mycoords = new rvec[natoms];
  }
  close_xtc(fit_file);
  delete [] mycoords;
  mycoords = NULL;
  cout << "done." << endl;

  // Open output files
  distances.open(d_filename.c_str());
  indices.open(i_filename.c_str());

  // Allocate vectors for storing the RMSDs for a structure
  fits.resize(ref_coords->size());
  permutation.resize(ref_coords->size());

  // Fix k if number of frames is too small
  if (ref_coords->size()-1 < k)
    k = ref_coords->size()-1;
  k1 = k + 1;
  keepers.resize(k1);

  // Timer for ETA
  time_t start = std::time(0);
  time_t last = start;

  // Compute fits
  for (int fit_frame = 0; fit_frame < fit_coords->size(); fit_frame++) {
    
    // Update user of progress
    if (std::time(0) - last > update_interval) {
      last = std::time(0);
      time_t eta = start + ((last-start) * fit_coords->size() / fit_frame);
      cout << "\rFrame: " << fit_frame << ", will finish " 
	   << string(std::ctime(&eta)).substr(0,20);
      cout.flush();
    }
    
    // Do Work
#pragma omp parallel for
    for (int ref_frame = 0; ref_frame < ref_coords->size(); ref_frame++) {
      do_fit(natoms,weights,
	     (*fit_coords)[fit_frame],
	     (*ref_coords)[ref_frame]);
      fits[ref_frame] = rmsdev(natoms,weights,
			       (*fit_coords)[fit_frame],
			       (*ref_coords)[ref_frame]) * 10.0;
    }

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
    distances.write((char*) &(keepers[1]), (sizeof(double)/sizeof(char)) * k);
    indices.write((char*) &(permutation[1]), (sizeof(int)/sizeof(char)) * k);

  }

  cout << endl << endl;

  // Clean coordinates
  for (vector<coord_array>::iterator itr = ref_coords->begin();
       itr != ref_coords->end(); itr++) delete [] (*itr);
  for (vector<coord_array>::iterator itr = fit_coords->begin();
       itr != fit_coords->end(); itr++) delete [] (*itr);
  delete ref_coords;
  delete fit_coords;
  delete [] weights;

  return 0;
}
