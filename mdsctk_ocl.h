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
// Copyright (c) 2012-2015, Joshua L. Phillips.
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
// in the README & COPYING files - if they are missing, get the
// official version at cnls.lanl.gov/~jphillips/.
// 
// To help us fund MDSCTK development, we humbly ask that you cite the
// papers on the package - you can find them in the top README file.
// 
// For more info, check our website at
// http://www.cs.mtsu.edu/~jphillips/software.html
// 
//

#ifndef MDSCTK_OCL_H
#define MDSCTK_OCL_H

#include <CL/cl.hpp>
#include <utility>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

struct OCLDevice {
  cl::Platform platform;
  cl::Context context;
  cl::Device device;
  cl::CommandQueue queue;
  size_t global_memory_size;
};


// Grab OpenCL device
bool getOCLDevice(OCLDevice &device, int ocl_device_id = 0);

// Get all OpenCL devices
std::vector<OCLDevice> getOCLDevices();

// Kernel building functions
cl::Kernel buildKernelFromString(const std::string source_code, std::string kernel_name, OCLDevice &ocl_device);
std::string testKernelFromString(const std::string source_code, std::string kernel_name, OCLDevice &ocl_device);
cl::Kernel buildKernel(std::string filename, std::string kernel_name, OCLDevice &ocl_device);
std::string testKernel(std::string filename, std::string kernel_name, OCLDevice &ocl_device);

// Run a kernel
bool enqueueKernel(OCLDevice &device, cl::Kernel &kernel, size_t size, cl::Event *event = NULL);

// Get runtime from event
double getExecutionTime(cl::Event &event);

// Auxiliary function to print the device list
void printOCLDeviceList(const std::vector<OCLDevice> &ocl_devices);

// MDSCTK-OCL Built-in Kernels

extern const char* euclidean_distance_KernelSource;

#endif
