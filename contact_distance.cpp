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
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <math.h>

// OpenMP
#include <omp.h>

// GROMACS
#include <gromacs/tpxio.h>
#include <gromacs/xtcio.h>
#include <gromacs/index.h>

using namespace std;
typedef real (*coord_array)[3];

int main(int argc, char* argv[]) {
  
  if (argc != 5 && argc != 4) {
    cerr << endl;
    cerr << "Usage: " << argv[0] << " [# threads] [topology file] [xtc file] <ndx file>" << endl;
    cerr << "   Computes the atomic contacts for structures in the given" << endl;
    cerr << "   xtc file. A topology PDB file and atom index file should" << endl;
    cerr << "   be provided for determining the atoms to compare." << endl;
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
  int nframes= 0;
  int update_interval = 1;
  const char* top_filename = argv[2];
  const char* ref_filename = argv[3];
  const char* ndx_filename = NULL;
  t_fileio *ref_file;
  rvec *mycoords = NULL;
  gmx_bool bOK = 1;
  double *contact = NULL;
  vector<coord_array> *ref_coords = NULL;
  real *weights = NULL;
  int        gnx1,gnx2;
  atom_id    *index1,*index2;
  char       *grpname1,*grpname2;
  ofstream   index;
  ofstream   data;

  // Remove C stdout (silly GROMACS warnings going every which stream!)
  int myout = dup(1);
  dup2(2,1);

  // EPS
  double eps = 1.0;
  do { eps /= 2.0; } while (1.0 + (eps / 2.0) != 1.0);
  eps = sqrt(eps);

  // Check for topology
  if (argc == 5)
    ndx_filename = argv[4];

  // Setup threads
  omp_set_num_threads(nthreads);

  // Get number of atoms and check xtc
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

  // Get atom selections
  cerr << "Please select two (non-overlapping) groups for contact profiling..." << endl;
  get_index(&top.atoms,ndx_filename,1,&gnx1,&index1,&grpname1);
  cerr << endl;
  get_index(&top.atoms,ndx_filename,1,&gnx2,&index2,&grpname2);
  cerr << endl;

  cerr << "Total grid size is " << gnx1 << " x " << gnx2 << " = " << (gnx1*gnx2) << endl;

  // Read coordinates and weight-center all structures
  cerr << "Reading reference coordinates from file: " << ref_filename << " ... ";
  ref_coords = new vector<coord_array>;
  ref_file = open_xtc(ref_filename,"r");
  mycoords = new rvec[natoms];
  while (read_next_xtc(ref_file, natoms, &step, &time, box, mycoords, &prec, &bOK)) {
    ref_coords->push_back(mycoords);
    mycoords = new rvec[natoms];
  }
  close_xtc(ref_file);
  delete [] mycoords;
  mycoords = NULL;
  nframes = ref_coords->size();
  cerr << "done." << endl;

  // Allocate vectors for storing the distances for a structure
  contact = new double[gnx1*gnx2];
  weights = new real[gnx1*gnx2];
  for (int x = 0; x < natoms; x++) weights[x] = top.atoms.atom[x].m;

#pragma omp parallel for
  for (int i = 0; i < gnx1; i++)
    for (int j = 0; j < gnx2; j++) {
      weights[(i*gnx2)+j] = log(top.atoms.atom[index1[i]].m * top.atoms.atom[index2[j]].m)*3.0;
    }

  // Get update frequency
  cerr.precision(8);
  cerr.setf(ios::fixed,ios::floatfield);
  update_interval = ceil(sqrt(ref_coords->size()));

  // Restore C stdout.
  dup2(myout,1);

  index.open("index.dat");
  data.open("data.dat");

  // Compute fits
  for (int frame = 0; frame < nframes; frame++) {

    // Update user of progress
    if (frame % update_interval == 0) {
      cerr << "\rWorking: " 
	   << ((double) (frame)) / 
	((double) (ref_coords->size())) * 100.0 << "%";
      cerr.flush();
    }
    
    // Do Work
#pragma omp parallel for
    for (int i = 0; i < gnx1*gnx2; i++)
      contact[i] = 0.0;

#pragma omp parallel for
    for (int i = 0; i < gnx1; i++) {
      int ii = index1[i];
      for (int j = 0; j < gnx2; j++) {
	int jj = index2[j];
	double d = 0.0;
	for (int k = 0; k < 3; k++)
	  d += (((*ref_coords)[frame][ii][k] - (*ref_coords)[frame][jj][k]) *
		((*ref_coords)[frame][ii][k] - (*ref_coords)[frame][jj][k]));
	d = sqrt(d) * 10.0;
	d = 1.0 / (1.0 + exp(d - weights[(i*gnx2)+j]));
	if (d > eps)
	  contact[(i*gnx2)+j] = d;
      } // j
    } // i


    double sum = 0.0;
#pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < gnx1*gnx2; i++)
      sum += contact[i];

    sum = 1.0; // No normalization...
    int total = 0;
#pragma omp parallel for reduction(+:total)
    for (int i = 0; i < gnx1*gnx2; i++)
      if (contact[i] > 0) {
	contact[i] /= sum;
	total++;
      }

    index.write((char*) &total, sizeof(int) / sizeof(char));
    for (int i = 0; i < gnx1*gnx2; i++)
      if (contact[i] > 0.0) {
	index.write((char*) &i, sizeof(int) / sizeof(char));
	data.write((char*) &contact[i], sizeof(double) / sizeof(char));
      }

    // cerr << frame << " " << total << endl;

  } // frame

  cerr << "\rWorking: " << 100.0 << "%" << endl;

  index.close();
  data.close();

  // Clean coordinates
  for (vector<coord_array>::iterator itr = ref_coords->begin();
       itr != ref_coords->end(); itr++) delete [] (*itr);
  delete ref_coords;

  delete [] contact;
  delete [] weights;
  
  return 0;
}
