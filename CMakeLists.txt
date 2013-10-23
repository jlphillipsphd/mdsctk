#
# 
#                This source code is part of
# 
#                        M D S C T K
# 
#       Molecular Dynamics Spectral Clustering ToolKit
# 
#                        VERSION 1.2.0
# Written by Joshua L. Phillips.
# Copyright (c) 2013, Joshua L. Phillips.
# check out http://github.com/douradopalmares/mdsctk/ for more information.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# If you want to redistribute modifications, please consider that
# derived work must not be called official MDSCTK. Details are found
# in the README & LICENSE files - if they are missing, get the
# official version at github.com/douradopalmares/mdsctk/.
# 
# To help us fund MDSCTK development, we humbly ask that you cite
# the papers on the package - you can find them in the top README file.
# 
# For more info, check our website at http://github.com/douradopalmares/mdsctk/
# 
#

cmake_minimum_required (VERSION 2.6)
project (MDSCTK)
# The version number.
set (MDSCTK_VERSION_MAJOR 1)
set (MDSCTK_VERSION_MINOR 2)

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (	
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
  @ONLY)
configure_file (	
  "${PROJECT_SOURCE_DIR}/config.r.in"
  "${PROJECT_BINARY_DIR}/config.r"
  @ONLY)
configure_file (	
  "${PROJECT_SOURCE_DIR}/MDSCTK.bash.in"
  "${PROJECT_BINARY_DIR}/MDSCTK.bash"
  @ONLY)

# add the binary tree to the search path for include files
# so that we will find config.h
include_directories("${PROJECT_BINARY_DIR}")

enable_language(CXX)

# Set to 'Debug' if no optimizations are wanted...
set(CMAKE_BUILD_TYPE Debug)

# Compiler flags
if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS_RELEASE
      "-O3 -Wl,-as-needed")  # Optimize
endif(CMAKE_COMPILER_IS_GNUCXX)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}/mdsctk" CACHE
    STRING "Installation prefix (installation will need write permissions here)" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

## Header files
include(CheckIncludeFileCXX)

## Oracle Berkeley DB
set(DB_CXX_INCLUDE_PATH ""
  CACHE STRING "Location of db_cxx.h")
unset(HAVE_DB_CXX_H CACHE)
if(NOT DB_CXX_INCLUDE_PATH STREQUAL "")
set(CMAKE_REQUIRED_FLAGS "-I${DB_CXX_INCLUDE_PATH}")
endif(NOT DB_CXX_INCLUDE_PATH STREQUAL "")
check_include_file_cxx(db_cxx.h      HAVE_DB_CXX_H)
if(NOT HAVE_DB_CXX_H)
  message(FATAL_ERROR "ERROR: unable to find Berekely DB C++ header file.")
endif(NOT HAVE_DB_CXX_H)

## OpenCL
set(OPENCL_INCLUDE_PATH ""
  CACHE STRING "Location of CL/cl.hpp")
unset(HAVE_OPENCL_H CACHE)
if(NOT OPENCL_INCLUDE_PATH STREQUAL "")
set(CMAKE_REQUIRED_FLAGS "-I${OPENCL_INCLUDE_PATH}")
endif(NOT OPENCL_INCLUDE_PATH STREQUAL "")
check_include_file_cxx(CL/cl.hpp      HAVE_OPENCL_H)
if(NOT HAVE_OPENCL_H)
  message(STATUS "WARNING: unable to find OpenCL header file.")
endif(NOT HAVE_OPENCL_H)

## Basic library checks
include(CheckLibraryExists)

## ARPACK
set(ARPACK_LIBRARY_PATH ""
  CACHE STRING "Location of libarpack.so")
unset(HAVE_LIBARPACK CACHE)
if(NOT ARPACK_LIBRARY_PATH STREQUAL "")
set(CMAKE_REQUIRED_FLAGS "-L${ARPACK_LIBRARY_PATH}")
endif(NOT ARPACK_LIBRARY_PATH STREQUAL "")
check_library_exists(arpack dsaupd_ "" HAVE_LIBARPACK)
if(HAVE_LIBARPACK)
  set(LIBS ${LIBS} arpack)
  if(NOT ARPACK_LIBRARY_PATH STREQUAL "")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${ARPACK_LIBRARY_PATH}")
  endif(NOT ARPACK_LIBRARY_PATH STREQUAL "")
else(HAVE_LIBARPACK)
  message(FATAL_ERROR "ERROR: unable to find ARPACK library.")
endif(HAVE_LIBARPACK)

## Oracle Berkeley DB
set(DB_CXX_LIBRARY_PATH ""
  CACHE STRING "Location of libdb_cxx.so")
if(NOT DB_CXX_LIBRARY_PATH STREQUAL "")
set(CMAKE_REQUIRED_FLAGS "-L${DB_CXX_LIBRARY_PATH}")
endif(NOT DB_CXX_LIBRARY_PATH STREQUAL "")
unset(HAVE_LIBDB_CXX CACHE)
check_library_exists(db_cxx open "" HAVE_LIBDB_CXX)
if(HAVE_LIBDB_CXX)
  set(LIBS ${LIBS} db_cxx)
  if(NOT DB_CXX_LIBRARY_PATH STREQUAL "")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${DB_CXX_LIBRARY_PATH}")
  endif(NOT DB_CXX_LIBRARY_PATH STREQUAL "")    
else(HAVE_LIBDB_CXX)
  message(FATAL_ERROR "ERROR: unable to find Berkeley DB C++ library.")
endif(HAVE_LIBDB_CXX)

## OpenCL
set(OPENCL_LIBRARY_PATH ""
  CACHE STRING "Location of libOpenCL.so")
unset(HAVE_LIBOPENCL CACHE)
if(NOT OPENCL_LIBRARY_PATH STREQUAL "")
set(CMAKE_REQUIRED_FLAGS "-L${OPENCL_LIBRARY_PATH}")
endif(NOT OPENCL_LIBRARY_PATH STREQUAL "")
check_library_exists(OpenCL clGetPlatformIDs "" HAVE_LIBOPENCL)
if(HAVE_LIBOPENCL)
  set(LIBS ${LIBS} OpenCL)
  if(NOT OPENCL_LIBRARY_PATH STREQUAL "")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${OPENCL_LIBRARY_PATH}")
  endif(NOT OPENCL_LIBRARY_PATH STREQUAL "")
else(HAVE_LIBOPENCL)
  message(STATUS "WARNING: unable to find OpenCL library (libOpenCL.so).")
endif(HAVE_LIBOPENCL)

# PKGCONFIG checks
find_package(PkgConfig REQUIRED)
include(FindPackageHandleStandardArgs)

# OpenMP
find_package(OpenMP REQUIRED)
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

# LAPACK
find_package(LAPACK REQUIRED)
include_directories(${LAPACK_INCLUDE_DIRS})
set(LIBS ${LIBS} ${LAPACK_LIBRARIES})

# BLAS
find_package(BLAS REQUIRED)
include_directories(${BLAS_INCLUDE_DIRS})
set(LIBS ${LIBS} ${BLAS_LIBRARIES})

# Boost
find_package(Boost 1.36.0 COMPONENTS program_options)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  set(LIBS ${LIBS} ${Boost_LIBRARIES})
endif(Boost_FOUND)

# GROMACS-libgmx
pkg_check_modules(PC_LIBGMX REQUIRED libgmx)
set(LIBGMX_DEFINITIONS ${PC_LIBGMX_CFLAGS_OTHER})
find_path(LIBGMX_INCLUDE_DIR gromacs/tpxio.h
  HINTS ${PC_LIBGMX_INCLUDEDIR} ${PC_LIBGMX_INCLUDE_DIRS}
  PATH_SUFFIXES gromacs)
find_library(LIBGMX_LIBRARY NAMES gmx libgmx 
  HINTS ${PC_LIBGMX_LIBDIR} ${PC_LIBGMX_LIBRARY_DIRS})
set(LIBGMX_INCLUDE_DIRS ${LIBGMX_INCLUDE_DIR})
set(LIBGMX_LIBRARIES ${LIBGMX_LIBRARY})

find_package_handle_standard_args(libgmx DEFAULT_MSG
  LIBGMX_LIBRARY LIBGMX_INCLUDE_DIR)
mark_as_advanced(LIBGMX_INCLUDE_DIR LIBGMX_LIBRARY)

include_directories(${LIBGMX_INCLUDE_DIRS})
set(LIBS ${LIBS} ${LIBGMX_LIBRARIES})

# GROMACS-libgmxana
pkg_check_modules(PC_LIBGMXANA REQUIRED libgmxana)
set(LIBGMXANA_DEFINITIONS ${PC_LIBGMXANA_CFLAGS_OTHER})
find_path(LIBGMXANA_INCLUDE_DIR gromacs/tpxio.h
  HINTS ${PC_LIBGMXANA_INCLUDEDIR} ${PC_LIBGMXANA_INCLUDE_DIRS}
  PATH_SUFFIXES gromacs)
find_library(LIBGMXANA_LIBRARY NAMES gmxana libgmxana 
  HINTS ${PC_LIBGMXANA_LIBDIR} ${PC_LIBGMXANA_LIBRARY_DIRS})
set(LIBGMXANA_INCLUDE_DIRS ${LIBGMXANA_INCLUDE_DIR})
set(LIBGMXANA_LIBRARIES ${LIBGMXANA_LIBRARY})

find_package_handle_standard_args(libgmxana DEFAULT_MSG
  LIBGMXANA_LIBRARY LIBGMXANA_INCLUDE_DIR)
mark_as_advanced(LIBGMXANA_INCLUDE_DIR LIBGMXANA_LIBRARY)

include_directories(${LIBGMXANA_INCLUDE_DIRS})
set(LIBS ${LIBS} ${LIBGMXANA_LIBRARIES})

# Framework library
add_library(mdsctk STATIC mdsctk.cpp)

# Build tools and scripts...
set(MDSCTK_TOOLS
  auto_decomp_sparse
  auto_decomp_sparse_nystrom
  auto_heir_decomp_sparse
  bb_xtc_to_phipsi
  check_xtc
  contact_profile
  decomp_dense
  decomp_sparse
  decomp_sparse_nystrom
  dijkstra
  knn_data
  knn_data_sparse
  knn_rms
  make_sysparse
  make_gesparse
  phipsi_to_sincos
  rms_test
  split_xtc)

foreach(TOOL ${MDSCTK_TOOLS})
  add_executable(${TOOL} ${TOOL}.cpp)
  target_link_libraries(${TOOL} ${LIBS} mdsctk)
endforeach(TOOL ${MDSCTK_TOOLS})

## OpenCL tools...
if(HAVE_OPENCL_H AND HAVE_LIBOPENCL)

  add_library(mdsctk_ocl STATIC mdsctk_ocl.cpp)

  set(MDSCTK_OCL_TOOLS
    knn_data_ocl)

  foreach(TOOL ${MDSCTK_OCL_TOOLS})
    add_executable(${TOOL} ${TOOL}.cpp)
    target_link_libraries(${TOOL} ${LIBS} mdsctk mdsctk_ocl)
  endforeach(TOOL ${MDSCTK_OCL_TOOLS})
  
else(HAVE_OPENCL_H AND HAVE_LIBOPENCL)
  message(STATUS "WARNING:")
  message(STATUS "WARNING: OpenCL accelerated tools will NOT be built.")
  message(STATUS "WARNING:")
endif(HAVE_OPENCL_H AND HAVE_LIBOPENCL)
