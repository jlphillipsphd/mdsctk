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
#include "mdsctk_ocl.h"

int main(int argc, char* argv[]) {

  const char* program_name = "knn_data_ocl";
  bool optsOK = true;
  gmx::initForCommandLine(&argc,&argv);
  copyright(program_name);
  cout << "   Computes the k nearest neighbors of all pairs of" << endl;
  cout << "   vectors in the given binary data files." << endl;
  cout << endl;
  cout << "   Use -h or --help to see the complete list of options." << endl;
  cout << endl;

  // Option vars...
  int k = 0;
  int vector_size = 0;
  string ref_filename;
  string fit_filename;
  string d_filename;
  string i_filename;
  bool sort;
  int ocl_device_id;

  // Declare the supported options.
  po::options_description cmdline_options;
  po::options_description program_options("Program options");
  program_options.add_options()
    ("help,h", "show this help message and exit")
    ("knn,k", po::value<int>(&k), "Input:  K-nearest neighbors (int)")
    ("sort,s",po::value<bool>(&sort)->default_value(true),"Input:  Find K-nn,false=full distance matix (bool)")
    ("vector-size,v", po::value<int>(&vector_size), "Input:  Data vector length (int)")
    ("reference-file,r", po::value<string>(&ref_filename)->default_value("reference.pts"), "Input:  Reference data file (string:filename)")
    ("fit-file,f", po::value<string>(&fit_filename), "Input:  Fitting data file (string:filename)")
    ("distance-file,d", po::value<string>(&d_filename)->default_value("distances.dat"), "Output: K-nn distances file (string:filename)")
    ("index-file,i", po::value<string>(&i_filename)->default_value("indices.dat"), "Output: K-nn indices file (string:filename)")    
    ;
  po::options_description ocl_options("OpenCL options");
  ocl_options.add_options()
    ("device-id,c", po::value<int>(&ocl_device_id)->default_value(0), "Selected OpenCL device number")
    ("list-devices,l", "List available OpenCL devices, then exit")
    ;
  cmdline_options.add(program_options).add(ocl_options);

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
  po::notify(vm);    

  if (vm.count("help")) {
    cout << "usage: " << program_name << " [options]" << endl;
    cout << cmdline_options << endl;
    return 1;
  }
  // Setup OpenCL
  OCLDevice ocl_device;
  if (!getOCLDevice(ocl_device, ocl_device_id))
    return -1;

  if (vm.count("list-devices")) {
    return 0;
  }
  if (!vm.count("knn") && sort) {
    cout << "ERROR: --knn not supplied." << endl;
    cout << endl;
    optsOK = false;
  }
  if (!vm.count("vector-size")) {
    cout << "ERROR: --vector-size not supplied." << endl;
    cout << endl;
    optsOK = false;
  }
  if (!vm.count("fit-file"))
    fit_filename = ref_filename;

  if (!optsOK) {
    return -1;
  }

  cout << "Running with the following options:" << endl;
  cout << "knn =            " << k << endl;
  cout << "sort =           " << sort << endl;
  cout << "vector-size =    " << vector_size << endl;
  cout << "reference-file = " << ref_filename << endl;
  cout << "fit-file =       " << fit_filename << endl;
  cout << "distance-file =  " << d_filename << endl;
  cout << "index-file =     " << i_filename << endl;
  cout << endl;

  // Local vars...
  vector<float*> *ref_coords = NULL;
  vector<float*> *fit_coords = NULL;
  int k1 = k + 1;
  int update_interval = 1;
  double *keepers = NULL;
  ofstream distances;
  ofstream indices;
  float time = 0.0;
  permutation< ::real> fits;

  const std::string kernel_source(euclidean_distance_KernelSource);
  cl::Kernel kernel = buildKernelFromString(kernel_source,
				  string("dist"),ocl_device);

  ref_coords = new vector<float*>;
  fit_coords = new vector<float*>;

  // Read coordinates
  cout << "Reading reference coordinates from file: " << ref_filename << " ... ";
  ifstream myfile;
  myfile.open(ref_filename.c_str());
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
  cout << "done." << endl;

  cout << "Reading fitting coordinates from file: " << fit_filename << " ... ";
  myfile.open(fit_filename.c_str());
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
  cout << "done." << endl;

  // Open output files
  distances.open(d_filename.c_str());
  indices.open(i_filename.c_str());

  // Allocate vectors for storing the RMSDs for a structure
  fits.data.resize(ref_coords->size());

  // Fix k if number of frames is too small
  if (ref_coords->size()-1 < k)
    k = ref_coords->size()-1;
  k1 = k + 1;
  keepers = new double[k1];

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

  // Timers...
  time_t start = std::time(0);
  time_t last = start;
  
  // Compute fits
  for (int fit_frame = 0; fit_frame < fit_coords->size(); fit_frame++) {
    
    // Update user of progress
    if (std::time(0) - last > update_interval) {
      last = std::time(0);
      time_t eta = start + ((last-start) * fit_coords->size() / fit_frame);
      cout << "\rFrame: " << fit_frame << ", will finish " 
	   << string(std::ctime(&eta)).substr(0,20);
      cout.flush();
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
				       &fits.data.at(0));
    
    // Sort
    if (sort) {
      fits.sort(k1);
      for (int x = 0; x < k1; x++)
	keepers[x] = (double) fits.data[fits.indices[x]];
      
      // Write out closest k RMSD alignment scores and indices
      distances.write((char*) &(keepers[1]), (sizeof(double)/sizeof(char)) * k);
      indices.write((char*) &(fits.indices[1]), (sizeof(int)/sizeof(char)) * k);
    }
    else {
      distances.write((char*) &(fits.data[0]), (sizeof(double)/sizeof(char)) * ref_coords->size());
      indices.write((char*) &(fits.indices[0]), (sizeof(int)/sizeof(char)) * ref_coords->size());
    }
    
  }

  cout << endl << endl;

  cout << "OpenCL Device Execution Time: " << time << endl;
  cout << endl;

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
