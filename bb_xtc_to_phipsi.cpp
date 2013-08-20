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

// C
#include <math.h>

// C++
#include <iostream>
#include <fstream>

// GROMACS
#include <gromacs/xtcio.h>

// Local
#include "config.h"

const double RAD2DEG = 180.0 / M_PI;
using namespace std;

void crossprod(float C[],
	       float x1, float y1, float z1,
	       float x2, float y2, float z2) {
  C[0] = ((y1 * z2) - (z1 * y2));
  C[1] = ((z1 * x2) - (x1 * z2));
  C[2] = ((x1 * y2) - (y1 * x2));
  return;
}

float torsion(float pos1[], float pos2[], float pos3[], float pos4[], bool degrees) {
  float L[3], Lnorm;
  float R[3], Rnorm;
  float S[3];
  float angle;

  crossprod(L,
	    (pos2[0] - pos1[0]), (pos2[1] - pos1[1]), (pos2[2] - pos1[2]),
	    (pos3[0] - pos2[0]), (pos3[1] - pos2[1]), (pos3[2] - pos2[2]));
  crossprod(R,
	    (pos4[0] - pos3[0]), (pos4[1] - pos3[1]), (pos4[2] - pos3[2]),
	    (pos2[0] - pos3[0]), (pos2[1] - pos3[1]), (pos2[2] - pos3[2]));

  Lnorm = sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]);
  Rnorm = sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2]);

  crossprod(S, L[0], L[1], L[2], R[0], R[1], R[2]);
  angle = (L[0]*R[0] + L[1]*R[1] + L[2]*R[2]) / (Lnorm * Rnorm);

  if (angle > 1.0) angle = 1.0;
  if (angle < -1.0) angle = -1.0;

  angle = acos( angle );
  if (degrees)
    angle = angle * RAD2DEG;

  if ((S[0] * (pos3[0]-pos2[0]) + 
       S[1] * (pos3[1]-pos2[1]) +
       S[2] * (pos3[2]-pos2[2])) < 0 )
    angle = -angle;

  return angle;
}

int main(int argc, char* argv[]) {

  if (argc < 2) {
    cerr << endl;
    cerr << "   MDSCTK " << MDSCTK_VERSION_MAJOR << "." << MDSCTK_VERSION_MINOR << endl;
    cerr << "   Copyright (C) 2013 Joshua L. Phillips" << endl;
    cerr << "   MDSCTK comes with ABSOLUTELY NO WARRANTY; see LICENSE for details." << endl;
    cerr << "   This is free software, and you are welcome to redistribute it" << endl;
    cerr << "   under certain conditions; see README.md for details." << endl;
    cerr << endl;
    cerr << "Usage: " << argv[0] << " [xtc file]..." << endl;
    cerr << "   Convert the provided xtc file to phipsi angles and" << endl;
    cerr << "   write the results to standard output." << endl;
    cerr << endl;
    return -1;
  }

  // Main variables
  int natoms;
  const char* xtc_file = argv[1];
  rvec* mycoords = NULL;
  double* mymat = NULL;
  t_fileio *myfile;
  ofstream output;

  // Used when reading XTC files...
  int step = 1;
  float time = 0.0;
  matrix box;
  float prec = 0.001;
  int frames = 0;
  gmx_bool bOK = 1;

  // Get number of atoms and allocate data structures
  myfile = open_xtc(xtc_file,"r");
  read_first_xtc(myfile,&natoms, &step, &time, box, &mycoords, &prec, &bOK);
  close_xtc(myfile);
  mycoords = new rvec[natoms];
  mymat = new double[(2*(natoms/3)-2)];
  output.open("phipsi.dat");

  // Convert coordinates
  for (int x = 1; x < argc; x++) {
    myfile = open_xtc(xtc_file,"r");
    while (read_next_xtc(myfile, natoms, &step, &time, box, mycoords, &prec, &bOK)) {
      int i_mat = 0;    
      for (int x = 0; x < natoms-3;) {
	mymat[i_mat++] = (double) torsion(mycoords[x],mycoords[x+1],mycoords[x+2],mycoords[x+3],false);
	x += 2;
	
	mymat[i_mat++] = (double) torsion(mycoords[x],mycoords[x+1],mycoords[x+2],mycoords[x+3],false);
	x += 1;
      }
      output.write((char*) mymat, (sizeof(double) / sizeof(char)) * (2*(natoms/3)-2));   
      frames++;
    }
    close_xtc(myfile);
  }
  output.close();

  cerr << "Wrote " << frames
       << " vectors of length " << (2*(natoms/3)-2)
       << " (" << (frames*(2*(natoms/3)-2)) << " total values)." << endl;

  // Clean up
  delete [] mycoords;
  delete [] mymat;

  return 0;
}
