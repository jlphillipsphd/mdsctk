//
// 
//                This source code is part of
// 
//                        M D S C T K
// 
//       Molecular Dynamics Spectral Clustering ToolKit
// 
//                        VERSION 1.2.5
//
// Written by Joshua L. Phillips.
// Copyright (c) 2012-2016, Joshua L. Phillips.
// Check out http://www.cs.mtsu.edu/~jphillips/software.html for more
// information.
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of the
// License, or (at your option) any later version.
// 
// If you want to redistribute modifications, please consider that
// derived work must not be called official MDSCTK. Details are found
// in the README & LICENSE files - if they are missing, get the
// official version at github.com/jlphillipsphd/mdsctk/.
// 
// To help us fund MDSCTK development, we humbly ask that you cite the
// papers on the package - you can find them in the top README file.
// 
// For more info, check our website at
// http://www.cs.mtsu.edu/~jphillips/software.html
// 
//

// Local
#include "config.h"
#include "mdsctk.h"

int main(int argc, char* argv[]) {
  
  const char* program_name = "ca_xtc_to_thetaphi";
  bool optsOK = true;
  gmx::initForCommandLine(&argc,&argv);
  copyright(program_name);
  cout << "   Convert the provided XTC file to theta and phi angles and" << endl;
  cout << "   write the results to the selected output file." << endl;
  cout << "   Note that this code anticipates a single protein" << endl;
  cout << "   chain, and only the CA atoms to be present in" << endl;
  cout << "   the XTC file." << endl;
  cout << endl;
  cout << "   Use -h or --help to see the complete list of options." << endl;
  cout << endl;
  
  // Option vars...
  string xtc_filename;
  string output_filename;
  
  // Declare the supported options.
  po::options_description cmdline_options;
  po::options_description program_options("Program options");
  program_options.add_options()
    ("help,h", "show this help message and exit")
    ("xtc-file,x", po::value<string>(&xtc_filename)->default_value("traj.xtc"), "Input:  Trajectory file (string:filename)")
    ("output-file,o", po::value<string>(&output_filename)->default_value("thetaphi.dat"), "Output: Theta-Phi angle data file (string:filename)")
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

  if (!optsOK) {
    return -1;
  }
  
  cout << "Running with the following options:" << endl;
  cout << "xtc-file =    " << xtc_filename << endl;
  cout << "output-file = " << output_filename << endl;
  cout << endl;
  
  // Main variables
  int natoms;
  rvec* mycoords = NULL;
  double* mymat = NULL;
  t_fileio *myfile;
  ofstream output;
  
  // Used when reading XTC files...
  long int step = 1;
  float time = 0.0;
  matrix box;
  float prec = 0.001;
  int frames = 0;
  gmx_bool bOK = 1;

  // Get number of atoms and allocate data structures
  myfile = open_xtc(xtc_filename.c_str(),"r");
  read_first_xtc(myfile,&natoms, &step, &time, box, &mycoords, &prec, &bOK);
  close_xtc(myfile);
  mycoords = new rvec[natoms];
  mymat = new double[(2*natoms)-5];
  output.open(output_filename.c_str());
  
  // Convert coordinates
  myfile = open_xtc(xtc_filename.c_str(),"r");
  while (read_next_xtc(myfile, natoms, &step, &time, box, mycoords, &prec, &bOK)) {
    int i_mat = 0;
    for (i_mat = 0; i_mat < natoms-2; i_mat++)
      mymat[i_mat] = (double) theta(mycoords[i_mat],
				    mycoords[i_mat+1],
				    mycoords[i_mat+2],false);
    for (i_mat = 0; i_mat < natoms-3; i_mat++) {
      mymat[i_mat+natoms-2] = (double) torsion(mycoords[i_mat],
					       mycoords[i_mat+1],
					       mycoords[i_mat+2],
					       mycoords[i_mat+3],false);
    }
    output.write((char*) mymat, (sizeof(double) / sizeof(char)) * (natoms*2-5));
    frames++;
  }
  close_xtc(myfile);
  output.close();

  cout << "Wrote " << frames
       << " vectors of length " << (2*natoms-5)
       << " (" << (frames*(2*natoms-5)) << " total values)." << endl;
  cout << endl;

  // Clean up
  delete [] mycoords;
  delete [] mymat;

  return 0;
}
