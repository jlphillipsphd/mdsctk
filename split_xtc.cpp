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

// C++
#include <iostream>

// GROMACS
#include <gromacs/xtcio.h>

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
    cerr << "Usage: " << argv[0] << " [xtc file] [n]" << endl;
    cerr << "   Splits the provided xtc file into two separate." << endl;
    cerr << "   xtc files (landmarks.xtc and remainder.xtc) where" << endl;
    cerr << "   landmarks.xtc will contain every n-th frame from" << endl;
    cerr << "   the input xtc file, while remainder.xtc will contain" << endl;
    cerr << "   all other frames from the input trajectory." << endl;
    cerr << endl;
    return -1;
  }

  // Main variables
  int natoms;
  int nframes = 0;
  int lframes = 0;
  int rframes = 0;
  const char* my_filename = argv[1];
  int n = atoi(argv[2]);
  t_fileio *my_file;
  t_fileio *landmark_file;
  t_fileio *remainder_file;
  rvec* mycoords = NULL;

  // Used when reading XTC files...
  int step = 1;
  float time = 0.0;
  float start_time = 0.0;
  matrix box;
  float prec = 0.001;
  gmx_bool bOK = 1;

  // Get number of atoms and allocate data structures
  my_file = open_xtc(my_filename,"r");
  read_first_xtc(my_file,&natoms, &step, &time, box, &mycoords, &prec, &bOK);
  close_xtc(my_file);
  my_file = open_xtc(my_filename,"r");
  mycoords = new rvec[natoms];

  landmark_file = open_xtc("landmarks.xtc","w");
  remainder_file = open_xtc("remainder.xtc","w");
  
  // Convert coordinates
  while (read_next_xtc(my_file, natoms, &step, &time, box, mycoords, &prec, &bOK)) {
    if (nframes % n == 0) {
      write_xtc(landmark_file,natoms,step,time,box,mycoords,prec);
      lframes++;
    }
    else {
      write_xtc(remainder_file,natoms,step,time,box,mycoords,prec);
      rframes++;
    }
    nframes++;
    if (nframes == 1) {
      start_time = time;
    }
  }
 
  // Clean up
  close_xtc(my_file);
  close_xtc(landmark_file);
  close_xtc(remainder_file);

  cout << "XTC Statistics" << endl;
  cout << "N_Frames: " << nframes << endl;
  cout << "N_Atoms: " << natoms << endl;
  cout << "Start_Time: " << start_time << endl;
  cout << "End_Time: " << time << endl;
  cout << "Precision: " << prec << endl;

  cout << endl;
  cout << "Frame distribution..." << endl;
  cout << "Landmark N_Frames: " << lframes << endl;
  cout << "Remainder N_Frames: " << rframes << endl;

  delete [] mycoords;

  return 0;
}
