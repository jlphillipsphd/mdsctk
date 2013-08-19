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

// Local
#include "config.h"

using namespace std;

#define BLKSIZE 65536

int main(int argc, char* argv[]) {

  if (argc != 1) {
    cerr << endl;
    cerr << "   MDSCTK " << MDSCTK_VERSION_MAJOR << "." << MDSCTK_VERSION_MINOR << endl;
    cerr << endl;
    cerr << "Usage: " << argv[0] << endl;
    cerr << "   Convert the (binary-double) angles from phipsi.dat" << endl;
    cerr << "   to sin-cos euclidean coordinates and write the results" << endl;
    cerr << "   to sincos.dat." << endl;
    cerr << endl;
    return -1;
  }

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
  myin.open("phipsi.dat");
  myout.open("sincos.dat");

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
  
  cerr << "Wrote " << (input_length/8) 
       << " sin-cos pairs (" << (input_length/4) << " total values)." << endl;

  delete [] data;
  delete [] result;
  return 0;
}
