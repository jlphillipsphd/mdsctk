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
#include <numeric>
#include <algorithm>
#include <strings.h>
#include <stdlib.h>
#include <math.h>

// OpenCL
#include "ocl.h"

// Local
#include "config.h"

using namespace std;

#define STRIFY(str) #str
#define STRINGIFY(str) STRIFY(str)

vector<float> fits;

bool compare(int left, int right) {
  return fits[left] < fits[right];
}

const char *KernelSource = "\n"			\
  "__kernel void dist(const int n,          \n"	\
  "                   const int d,          \n" \
  "                __global const float *A, \n"	\
  "                __global const float *B, \n" \
  "                __global float *C) {     \n" \
  "                                         \n" \
  "  int i = get_global_id(0);              \n" \
  "  int j;                                 \n" \
  "  float x = 0.0;                         \n" \
  "                                         \n" \
  "  if (i < n) {                           \n" \
  "  for (j = 0; j < d; j++) {              \n" \
  "    x += (A[(i*d)+j] - B[j]) * (A[(i*d)+j] - B[j]);\n" \
  "  }                                      \n" \
  "  //printf(\"%d %f %f %f %f %f\\n\",i,A[(i*d)],A[(i*d)+1],B[0],B[1],sqrt(x));\n" \
  "  C[i] = sqrt(x);                        \n" \
  "  }                                      \n" \
  "}                                        \n" \
  "\n";

int main(int argc, char* argv[]) {

    cout << endl;
    cout << "   MDSCTK " << MDSCTK_VERSION_MAJOR << "." << MDSCTK_VERSION_MINOR << endl;
    cout << endl;
  
  if (argc != 5) {
    cerr << "Usage: " << argv[0] << " [k] [vector size] [reference data file] [fitting data file]" << endl;
    cerr << "   Computes the k nearest neighbors of all pairs of" << endl;
    cerr << "   vectors in the given binary data files." << endl;
    cerr << endl;
    return -1;
  }

  // Local vars
  int vector_size = 0;
  vector<float*> *ref_coords = NULL;
  vector<float*> *fit_coords = NULL;
  int k = atoi(argv[1]);
  int k1 = k + 1;
  vector_size = atoi(argv[2]);
  int update_interval = 1;
  const char* ref_file = argv[3];
  const char* fit_file = argv[4];
  double *keepers = NULL;
  ofstream distances;
  ofstream indices;
  float time = 0.0;
  vector<int> permutation;
  const std::string kernel_source(KernelSource);

  // Setup OpenCL
  OCLDevice ocl_device;
  if (!getOCLDevice(ocl_device))
    return -1;
  cout << "Using kernel file: " << string(STRINGIFY(MDSCTK_KERNELS))
    +string("/dist.cl") << endl;
  cl::Kernel kernel = buildKernelFromString(kernel_source,
				  string("dist"),ocl_device);
  // cl::Kernel kernel = buildKernel(string(STRINGIFY(MDSCTK_KERNELS))+string("/dist.cl"),
  // 				  string("dist"),ocl_device);

  ref_coords = new vector<float*>;
  fit_coords = new vector<float*>;

  // Read coordinates
  cerr << "Reading reference coordinates from file: " << ref_file << " ... ";
  ifstream myfile;
  myfile.open(ref_file);
  double* myread = new double[vector_size];
  float* mycoords = new float[vector_size];
  myfile.read((char*) myread, sizeof(double) * vector_size);
  while (!myfile.eof()) {
    for (int x = 0; x < vector_size; x++)
      mycoords[x] = (float) myread[x];
    ref_coords->push_back(mycoords);
    mycoords = new float[vector_size];
    myfile.read((char*) myread, sizeof(double) * vector_size);
  }
  myfile.close();
  cerr << "done." << endl;

  cerr << "Reading fitting coordinates from file: " << fit_file << " ... ";
  myfile.open(fit_file);
  myfile.read((char*) myread, sizeof(double) * vector_size);
  while (!myfile.eof()) {
    for (int x = 0; x < vector_size; x++)
      mycoords[x] = (float) myread[x];
    fit_coords->push_back(mycoords);
    mycoords = new float[vector_size];
    myfile.read((char*) myread, sizeof(double) * vector_size);
  }
  myfile.close();
  delete [] mycoords;
  mycoords = NULL;
  cerr << "done." << endl;

  // Open output files
  distances.open("distances.dat");
  indices.open("indices.dat");

  // Allocate vectors for storing the RMSDs for a structure
  fits.resize(ref_coords->size());
  permutation.resize(ref_coords->size());

  // Fix k if number of frames is too small
  if (ref_coords->size()-1 < k)
    k = ref_coords->size()-1;
  k1 = k + 1;
  keepers = new double[k1];

  // Get update frequency
  // int update_interval = (int) floor(sqrt((float) coords.size()));
  cerr.precision(8);
  cerr.setf(ios::fixed,ios::floatfield);
  update_interval = ceil(sqrt(fit_coords->size()));

  // Setup OCL
  cl::Buffer bufferA = cl::Buffer(ocl_device.context, CL_MEM_READ_ONLY,
				  ref_coords->size() * vector_size * sizeof(float));
  cl::Buffer bufferB = cl::Buffer(ocl_device.context, CL_MEM_READ_ONLY,
				  vector_size * sizeof(float));
  cl::Buffer bufferC = cl::Buffer(ocl_device.context, CL_MEM_WRITE_ONLY,
				  ref_coords->size() * sizeof(float));
  kernel.setArg(0, (int) ref_coords->size());
  kernel.setArg(1, (int) vector_size);
  kernel.setArg(2, bufferA);
  kernel.setArg(3, bufferB);
  kernel.setArg(4, bufferC);
  
  // Load reference frames...
  for (int ref_frame = 0; ref_frame < ref_coords->size(); ref_frame++)
    ocl_device.queue.enqueueWriteBuffer(bufferA, CL_TRUE,
					ref_frame * vector_size * sizeof(float),
					vector_size * sizeof(float),
					(*ref_coords)[ref_frame]);

  // Compute fits
  for (int fit_frame = 0; fit_frame < fit_coords->size(); fit_frame++) {
    
    // Update user of progress
    if (fit_frame % update_interval == 0) {
      cerr << "\rWorking: " << (((float) fit_frame) / ((float) fit_coords->size())) * 100.0 << "%";
      cerr.flush();
    }

    // Do Work
    cl::Event event;
    ocl_device.queue.enqueueWriteBuffer(bufferB, CL_TRUE,
					0,
					vector_size * sizeof(float),
					(*fit_coords)[fit_frame]);
    if (!enqueueKernel(ocl_device,kernel,ref_coords->size(),&event))
      return -1;
    ocl_device.queue.finish();
    time += getExecutionTime(event);
    ocl_device.queue.enqueueReadBuffer(bufferC, CL_TRUE,
				       0,
				       ref_coords->size() * sizeof(float),
				       &fits.at(0));

    // Sort
    for (int x = 0; x < permutation.size(); x++)
      permutation[x] = x;
    sort(permutation.begin(), permutation.end(), compare);
    for (int x = 0; x < k1; x++)
      keepers[x] = (double) fits[permutation[x]];
    
    // Write out closest k RMSD alignment scores and indices
    distances.write((char*) &(keepers[1]), (sizeof(double) / sizeof(char)) * k);
    indices.write((char*) &(permutation[1]), (sizeof(int) / sizeof(char)) * k);
  }

  cerr << "\rWorking: " << 100.0 << "%" << endl;

  cerr << "OpenCL Device Execution Time: " << time << endl;
  cerr << endl;

  // Clean coordinates
  for (vector<float*>::iterator itr = ref_coords->begin();
       itr != ref_coords->end(); itr++) delete [] (*itr);
  for (vector<float*>::iterator itr = fit_coords->begin();
       itr != fit_coords->end(); itr++) delete [] (*itr);
  delete ref_coords;
  delete fit_coords;
  delete [] keepers;
  delete [] myread;

  return 0;
}
