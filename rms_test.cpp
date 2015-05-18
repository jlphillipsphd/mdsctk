//
// 
//                This source code is part of
// 
//                        M D S C T K
// 
//       Molecular Dynamics Spectral Clustering ToolKit
// 
//                        VERSION 1.2.2
// Written by Joshua L. Phillips.
// Copyright (c) 2012-2014, Joshua L. Phillips.
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

// Local
#include "config.h"
#include "mdsctk.h"

int main(int argc, char* argv[]) {

  const char* program_name = "rms_test";
  bool optsOK = true;
  gmx::initForCommandLine(&argc,&argv);
  copyright(program_name);
  cout << "   Computes the RMSD of all structures for the given" << endl;
  cout << "   in the xtc file. The topology file provided will be" << endl;
  cout << "   used as the reference structure." << endl;
  cout << endl;
  cout << "   Use -h or --help to see the complete list of options." << endl;
  cout << endl;

  // Option vars...
  string top_filename;
  string xtc_filename;
  string output_filename;

  // Declare the supported options.
  po::options_description cmdline_options;
  po::options_description program_options("Program options");
  program_options.add_options()
    ("help,h", "show this help message and exit")
    ("topology-file,p", po::value<string>(&top_filename)->default_value("topology.pdb"), "Input:  Topology file [.pdb,.gro,.tpr] (string:filename)")
    ("xtc-file,x", po::value<string>(&xtc_filename)->default_value("traj.xtc"), "Input:  Trajectory file (string:filename)")
    ("output-file,o", po::value<string>(&output_filename)->default_value("rms.dat"), "Output: RMS data file (string:filename)")    
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
  cout << "topology-file = " << top_filename << endl;
  cout << "xtc-file =      " << xtc_filename << endl;
  cout << "output-file =   " << output_filename << endl;
  cout << endl;
  
  coord_array fit_coords = NULL;
  ofstream output;

  // Get number of atoms and initialize weights
  cout << "Reading topology coordinates from file: " << top_filename << " ... ";
  TOP_file ref_file(top_filename);
  cout << "done." << endl;
  cout << "Reading fitting coordinates from file: " << xtc_filename << " ... ";
  XTC_file fit_file(xtc_filename);

  if (fit_file.get_natoms() != ref_file.get_natoms()) {
    cout << "*** ERROR ***" << endl;
    cout << "Number of atoms in topology file ("
	 << ref_file.get_natoms() << ") "
	 << "does not match the number of atoms "
	 << "in the XTC file (" << fit_file.get_natoms() << ")."
	 << endl;
    exit(4);
  }

  // Output
  output.open(output_filename.c_str());

  // Read coordinates and weight-center all structures
  fit_coords = fit_file.get_next_frame_ptr();
  while (fit_coords) {
    ref_file.center(fit_coords);
    output << fit_file.get_time() << "\t" 
	   << ref_file.rmsd(ref_file.get_frame_ptr(),fit_coords) << endl;
    fit_coords = fit_file.get_next_frame_ptr();
  }
  output.close();
  cout << "done." << endl;

  return 0;
}
