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
// official version at github.com/douradopalmares/mdsctk/.
// 
// To help us fund MDSCTK development, we humbly ask that you cite the
// papers on the package - you can find them in the top README file.
// 
// For more info, check our website at
// http://www.cs.mtsu.edu/~jphillips/software.html
// 
//

#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>
#include <utility>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "mdsctk_ocl.h"

std::vector<OCLDevice> getOCLDevices() {
  std::cerr << "Initializing OpenCL..." << std::endl;
  std::cerr << std::endl;
  std::vector<OCLDevice> ocl_devices;
  try { 
    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);
    for (int p = 0 ; p < platforms.size(); p++) {
      cl_context_properties cps[3] = { 
	CL_CONTEXT_PLATFORM, 
	(cl_context_properties)(platforms[p])(), 
	0 
      };
      cl::Context context( CL_DEVICE_TYPE_DEFAULT, cps);
      std::vector<cl::Device> devices = context.getInfo<CL_CONTEXT_DEVICES>();
      for (int d = 0; d < devices.size(); d++) {
	OCLDevice ocl_device;
	// Create a command queue and use the first device
	ocl_device.platform = platforms[p];
	ocl_device.device = devices[d];
	ocl_device.queue = cl::CommandQueue(context, devices[d], CL_QUEUE_PROFILING_ENABLE);
	ocl_device.context = context;
	ocl_device.global_memory_size = ocl_device.device.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>();
	ocl_devices.push_back(ocl_device);
      }
    }
  } catch(cl::Error error) {
    std::cerr << "OpenCL Exception: " << error.what() << "(" << error.err() << ")" << std::endl;
  }

  if (ocl_devices.size() == 0) {
    std::cerr << "Could NOT locate an OpenCL platform..." << std::endl;
    std::cerr << std::endl;
    std::cerr << "Please make sure that you have installed" << std::endl;
    std::cerr << "the .icd files for your platforms in" << std::endl;
    std::cerr << "/etc/OpenCL/vendors or have set the" << std::endl;
    std::cerr << "OCL_ICD_VENDORS path appropriately..." << std::endl;
    std::cerr << std::endl;
  }

  return ocl_devices;
}

void printOCLDeviceList(const std::vector<OCLDevice> &ocl_devices) {
  std::cout << "OpenCL Device List" << std::endl;
  std::cout << std::endl;
  for (int d = 0; d < ocl_devices.size(); d++) {
    std::cout << "[" << d << "]" << std::endl;
    std::cout << "     Vendor: " << ocl_devices[d].platform.getInfo<CL_PLATFORM_VENDOR>() << std::endl;
    std::cout << "     Name: " << ocl_devices[d].device.getInfo<CL_DEVICE_NAME>() << std::endl;
    std::cout << "     Global Memory: " << (ocl_devices[d].global_memory_size / 1024 / 1024) << "MB" << std::endl;
  }
  std::cout << std::endl;
}

bool getOCLDevice(OCLDevice &device, int ocl_device_id) {
  std::vector<OCLDevice> ocl_devices = getOCLDevices();
  
  if (ocl_devices.size() == 0)
    return false;

  printOCLDeviceList(ocl_devices);
  
  std::cout << "Selected OpenCL Device: [" << ocl_device_id << "]" << std::endl;
  std::cout << std::endl;

  if (ocl_device_id >= ocl_devices.size()) {
    std::cerr << "Selected Illegal Device [" << ocl_device_id << "]" << std::endl; 
    std::cerr << std::endl;
    return false;
  }

  device = ocl_devices[ocl_device_id];
  return true;
}

cl::Kernel buildKernelFromString(const std::string source_code, std::string kernel_name, OCLDevice &ocl_device) {
  cl::Kernel kernel;
  try { 
    // Read source string
    cl::Program::Sources source(1, std::make_pair(source_code.c_str(), source_code.length()+1));
    
    // Make program of the source code in the context
    cl::Program program = cl::Program(ocl_device.context, source);
    std::vector<cl::Device> devices = ocl_device.context.getInfo<CL_CONTEXT_DEVICES>();
    
    // Build program for these specific devices
    program.build(devices);
    
    // Make kernel
    kernel = cl::Kernel(program, kernel_name.c_str());
  } catch(cl::Error error) {
    std::cerr << "OpenCL Exception: " << error.what() << "(" << error.err() << ")" << std::endl;
  }
  return kernel;
}

std::string testKernelFromString(const std::string source_code, std::string kernel_name, OCLDevice &ocl_device) {
  cl::Kernel kernel;
  cl::Program program;
  std::string log = "SUCCESS!!!";
  try { 
    // Read source string
    cl::Program::Sources source(1, std::make_pair(source_code.c_str(), source_code.length()+1));

    // Make program of the source code in the context
    program = cl::Program(ocl_device.context, source);
    std::vector<cl::Device> devices = ocl_device.context.getInfo<CL_CONTEXT_DEVICES>();

    // Build program for these specific devices
    program.build(devices);

    // Make kernel
    kernel = cl::Kernel(program, kernel_name.c_str());
    
  } catch(cl::Error error) {
    std::cerr << "OpenCL Exception: " << error.what() << "(" << error.err() << ")" << std::endl;
    log = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(ocl_device.device);
  }
  return log;
}

cl::Kernel buildKernel(std::string filename, std::string kernel_name, OCLDevice &ocl_device) {
  cl::Kernel kernel;
  try { 
    // Read source file
    std::ifstream sourceFile(filename.c_str());
    std::string sourceCode(std::istreambuf_iterator<char>(sourceFile),
			   (std::istreambuf_iterator<char>()));
    cl::Program::Sources source(1, std::make_pair(sourceCode.c_str(), sourceCode.length()+1));
    
    // Make program of the source code in the context
    cl::Program program = cl::Program(ocl_device.context, source);
    std::vector<cl::Device> devices = ocl_device.context.getInfo<CL_CONTEXT_DEVICES>();
    
    // Build program for these specific devices
    program.build(devices);
    
    // Make kernel
    kernel = cl::Kernel(program, kernel_name.c_str());
  } catch(cl::Error error) {
    std::cerr << "OpenCL Exception: " << error.what() << "(" << error.err() << ")" << std::endl;
  }
  return kernel;
}


std::string testKernel(std::string filename, std::string kernel_name, OCLDevice &ocl_device) {
  cl::Kernel kernel;
  cl::Program program;
  std::string log = "SUCCESS!!!";
  try { 
    // Read source file
    std::ifstream sourceFile(filename.c_str());
    std::string sourceCode(std::istreambuf_iterator<char>(sourceFile),
			   (std::istreambuf_iterator<char>()));
    cl::Program::Sources source(1, std::make_pair(sourceCode.c_str(), sourceCode.length()+1));

    // Make program of the source code in the context
    program = cl::Program(ocl_device.context, source);
    std::vector<cl::Device> devices = ocl_device.context.getInfo<CL_CONTEXT_DEVICES>();

    // Build program for these specific devices
    program.build(devices);

    // Make kernel
    kernel = cl::Kernel(program, kernel_name.c_str());
    
  } catch(cl::Error error) {
    std::cerr << "OpenCL Exception: " << error.what() << "(" << error.err() << ")" << std::endl;
    log = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(ocl_device.device);
  }
  return log;
}

bool enqueueKernel(OCLDevice &device, cl::Kernel &kernel, size_t size, cl::Event *event) {
  try { 
    // Run the kernel on specific ND range
    cl::NDRange local(kernel.getWorkGroupInfo<CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE>(device.device));
    cl::NDRange global(size + (local[0] - (size % local[0])));
    if (event)
      device.queue.enqueueNDRangeKernel(kernel, cl::NullRange, global, local, NULL, event);
    else
      device.queue.enqueueNDRangeKernel(kernel, cl::NullRange, global, local);
  } catch(cl::Error error) {
    std::cerr << "OpenCL Exception: " << error.what() << "(" << error.err() << ")" << std::endl;
    return false;
  } 
  return true;
}

double getExecutionTime(cl::Event &event) {
  try { 
    size_t n_begin, n_end;
    event.wait();
    n_begin = event.getProfilingInfo<CL_PROFILING_COMMAND_START>();
    n_end = event.getProfilingInfo<CL_PROFILING_COMMAND_END>();
    return 1.0E-9 * (n_end - n_begin);
  } catch(cl::Error error) {
    std::cerr << "OpenCL Exception: " << error.what() << "(" << error.err() << ")" << std::endl;
    return -0.0;
  } 
}

const char *euclidean_distance_KernelSource = "\n" \
  "#define AX(i,j,d) A[(i*d)+j]             \n" \
  "                                         \n" \
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
  "   if (i < n) {                          \n" \
  "      for (j = 0; j < d; j++) {          \n" \
  "        x += (AX(i,j,d) - B[j]) * (AX(i,j,d) - B[j]);\n" \
  "      }                                  \n" \
  "      C[i] = sqrt(x);                    \n" \
  "   }                                     \n" \
  "}                                        \n" \
  "\n";
