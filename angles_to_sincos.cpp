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
// official version at github.com/douradopalmares/mdsctk/.
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

  const char* program_name = "angles_to_sincos";
  bool optsOK = true;
  gmx::initForCommandLine(&argc,&argv);
  copyright(program_name);
  cout << "   Convert the (binary-double) angles from the input file" << endl;
  cout << "   to sin-cos euclidean coordinates and write the results" << endl;
  cout << "   to the provided output file." << endl;
  cout << endl;
  cout << "   Use -h or --help to see the complete list of options." << endl;
  cout << endl;

  // Option vars...
  string input_filename;
  string output_filename;

  // Declare the supported options.
  po::options_description cmdline_options;
  po::options_description program_options("Program options");
  program_options.add_options()
    ("help,h", "show this help message and exit")
    ("input-file,i", po::value<string>(&input_filename)->default_value("phipsi.dat"), "Input:  Phi-psi angle data file (string:filename)")
    ("output-file,o", po::value<string>(&output_filename)->default_value("sincos.dat"), "Output: Projected angle data file (string:filename)")    
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
  cout << "input-file  = " << input_filename << endl;
  cout << "output-file = " << output_filename << endl;
  cout << endl;

  int input_length = 0;
  int char_block_size = BLKSIZE / sizeof(char);
  int double_block_size = BLKSIZE / sizeof(double);
  double *data = new double[double_block_size];
  double *result = new double[double_block_size*2];
  int num_blocks = 0;
  int char_extra = 0;
  int double_extra = 0;

  ifstream myin;
  ofstream myout;
  myin.open(input_filename.c_str());
  myout.open(output_filename.c_str());

  myin.seekg(0, ios::end);
  input_length = myin.tellg();
  myin.seekg(0, ios::beg);

  num_blocks = input_length / char_block_size;
  char_extra = input_length % char_block_size;
  double_extra = (input_length % char_block_size) / sizeof(double);

  for (int x = 0; x < num_blocks; x++) {
    myin.read((char*) data, char_block_size);
    for (int y = 0; y < double_block_size; y++) {
      result[2*y] = sin(data[y]);
      result[(2*y)+1] = cos(data[y]);
    } 
    myout.write((char*) result, char_block_size*2);
  }

  myin.read((char*) data, char_extra);
  for (int y = 0; y < double_extra; y++) {
    result[2*y] = sin(data[y]);
    result[(2*y)+1] = cos(data[y]);
  } 
  myout.write((char*) result, char_extra * 2);
  
  cout << "Wrote " << (input_length/8) 
       << " sin-cos pairs (" << (input_length/4) << " total values)." << endl;
  cout << endl;

  delete [] data;
  delete [] result;
  return 0;
}
