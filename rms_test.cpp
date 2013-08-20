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
// C
#include <string.h>
#include <stdlib.h>
#include <math.h>
// C++
#include <iostream>
#include <fstream>
#include <vector>

// GROMACS
#include <gromacs/tpxio.h>
#include <gromacs/xtcio.h>
#include <gromacs/do_fit.h>

// Local
#include "config.h"

using namespace std;

int main(int argc, char* argv[]) {
  
  if (argc != 3) {
    cerr << endl;
    cerr << "   MDSCTK " << MDSCTK_VERSION_MAJOR << "." << MDSCTK_VERSION_MINOR << endl;
    cerr << "   Copyright (C) 2013 Joshua L. Phillips" << endl;
    cerr << "   MDSCTK comes with ABSOLUTELY NO WARRANTY; see LICENSE for details." << endl;
    cerr << "   This is free software, and you are welcome to redistribute it" << endl;
    cerr << "   under certain conditions; see README.md for details." << endl;
    cerr << endl;
    cerr << "Usage: " << argv[0] << " [reference structure] [fitting xtc file]" << endl;
    cerr << "   Computes the RMSD of all structures for the given" << endl;
    cerr << "   in the xtc file. A template structure should be" << endl;
    cerr << "   provided as the reference structure." << endl;
    cerr << endl;
    return -1;
  }

  int ref_natoms = -1;
  int natoms = 0;
  rvec* ref_coords = NULL;
  rvec* fit_coords = NULL;
  real *weights = NULL;
  gmx_bool bOK = 1;

  const char *ref_filename = argv[1];
  const char *fit_filename = argv[2];
  t_fileio *fit_file;

  // Used when reading XTC files...
  int step = 1;
  float time = 0.0;
  matrix box;
  float prec = 0.001;
  char buf[256];
  t_topology top;
  int ePBC;

  // Get number of atoms and initialize weights
  cerr << "Reading reference coordinates from " << ref_filename << " ... ";
  read_tps_conf(ref_filename, buf, &top, &ePBC, &ref_coords,
		NULL, box, TRUE);
  ref_natoms = top.atoms.nr;
  cerr << "done." << endl;

  fit_file = open_xtc(fit_filename,"r");
  read_first_xtc(fit_file,&natoms, &step, &time, box, &fit_coords, &prec, &bOK);
  close_xtc(fit_file);
  if (natoms != ref_natoms) {
    cerr << "*** ERROR ***" << endl;
    cerr << "Number of atoms in topology file ("
	 << ref_natoms << ") "
	 << "does not match the number of atoms "
	 << "in the XTC file (" << natoms << ")."
	 << endl;
    exit(4);
  }

  fit_coords = new rvec[natoms];
  weights = new real[natoms];
  for (int x = 0; x < natoms; x++) {
    weights[x] = top.atoms.atom[x].m;
  }
  reset_x(natoms,NULL,natoms,NULL,ref_coords,weights);

  // Read coordinates and weight-center all structures
  cerr << "Reading fitting coordinates from file: " << fit_filename << " ... ";
  fit_file = open_xtc(fit_filename,"r");
  while (read_next_xtc(fit_file, natoms, &step, &time, box, fit_coords, &prec, &bOK)) {
    reset_x(natoms,NULL,natoms,NULL,fit_coords,weights);
    do_fit(natoms,weights,ref_coords,fit_coords);
    cout << time << "\t" << (rmsdev(natoms,weights,ref_coords,fit_coords) * 10.0) << endl;
  }
  close_xtc(fit_file);
  cerr << "done." << endl;

  // Clean coordinates
  delete [] ref_coords;
  delete [] fit_coords;
  delete [] weights;
  
  return 0;
}
