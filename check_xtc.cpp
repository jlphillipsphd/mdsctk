//
// 
//                This source code is part of
// 
//                        M D S C T K
// 
//       Molecular Dynamics Spectral Clustering ToolKit
// 
//                        VERSION 1.2.5
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

  const char* program_name = "check_xtc";
  bool optsOK = true;
  gmx::initForCommandLine(&argc,&argv);
  copyright(program_name);
  cout << "   Report stats on the provided xtc file." << endl;
  cout << endl;
  cout << "   Use -h or --help to see the complete list of options." << endl;
  cout << endl;

  // Option vars...
  string xtc_filename;

  // Declare the supported options.
  po::options_description cmdline_options;
  po::options_description program_options("Program options");
  program_options.add_options()
    ("help,h", "show this help message and exit")
    ("xtc-file,x", po::value<string>(&xtc_filename)->default_value("traj.xtc"), "Input:  Trajectory file (string:filename)")
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
  cout << endl;

  // Main variables
  int natoms;
  int nframes = 0;
  t_fileio *my_file;
  rvec* mycoords = NULL;

  // Used when reading XTC files...
  long int step = 1;
  float time = 0.0;
  float start_time = 0.0;
  matrix box;
  float prec = 0.001;
  gmx_bool bOK = 1;

  // Get number of atoms and allocate data structures
  my_file = open_xtc(xtc_filename.c_str(),"r");
  read_first_xtc(my_file,&natoms, &step, &time, box, &mycoords, &prec, &bOK);
  close_xtc(my_file);
  my_file = open_xtc(xtc_filename.c_str(),"r");
  mycoords = new rvec[natoms];
  
  // Convert coordinates
  while (read_next_xtc(my_file, natoms, &step, &time, box, mycoords, &prec, &bOK)) {
    nframes++;
    if (nframes == 1) {
      start_time = time;
    }
  }
 
  // Clean up
  close_xtc(my_file);

  cout << "XTC Statistics - " << xtc_filename << endl;
  cout << "Number of frames: " << nframes << endl;
  cout << "Nunber of atoms:  " << natoms << endl;
  cout << "Start time:       " << start_time << endl;
  cout << "End time:         " << time << endl;
  cout << "Precision:        " << prec << endl;
  cout << endl;

  delete [] mycoords;

  return 0;
}
