#
# 
#                This source code is part of
# 
#                        M D S C T K
# 
#       Molecular Dynamics Spectral Clustering ToolKit
# 
#                        VERSION 1.1.1
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
set (MDSCTK_VERSION_MINOR 1)

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (	
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
  )					
configure_file (	
  "${PROJECT_SOURCE_DIR}/config.r.in"
  "${PROJECT_BINARY_DIR}/config.r"
  )					

# add the binary tree to the search path for include files
# so that we will find config.h
include_directories("${PROJECT_BINARY_DIR}")

enable_language(CXX)

# Set to 'Debug' if no optimizations are wanted...
set(CMAKE_BUILD_TYPE Debug)

# Compiler flags
if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS_RELEASE
      "-O3 -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF -Wl,-as-needed")  # Optimize
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-s")  # Strip binary
endif(CMAKE_COMPILER_IS_GNUCXX)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "/usr/local/mdsctk" CACHE
    STRING "Installation prefix (installation will need write permissions here)" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

## Header files
include(CheckIncludeFileCXX)

## Oracle Berkeley DB
set(DB_CXX_INCLUDE_PATH "/usr/include"
  CACHE STRING "Location of db_cxx.h")
unset(HAVE_DB_CXX_H CACHE)
set(CMAKE_REQUIRED_FLAGS "-I${DB_CXX_INCLUDE_PATH}")
check_include_file_cxx(db_cxx.h      HAVE_DB_CXX_H)
if(NOT HAVE_DB_CXX_H)
  message(FATAL_ERROR "ERROR: unable to find Berekely DB C++ header file.")
endif(NOT HAVE_DB_CXX_H)

## Boost
set(BOOST_INCLUDE_PATH "/usr/include"
  CACHE STRING "Location of boost/graph/dijkstra_shortest_paths.hpp")
unset(HAVE_BOOST_DSP_H CACHE)
set(CMAKE_REQUIRED_FLAGS "-I${BOOST_INCLUDE_PATH}")
check_include_file_cxx(boost/graph/dijkstra_shortest_paths.hpp
  HAVE_BOOST_DSP_H)
if(NOT HAVE_BOOST_DSP_H)
  message(FATAL_ERROR "ERROR: unable to find BOOST(C++) graph library header.")
endif(NOT HAVE_BOOST_DSP_H)

## OpenCL
set(OPENCL_INCLUDE_PATH "/usr/include"
  CACHE STRING "Location of CL/cl.hpp")
unset(HAVE_OPENCL_H CACHE)
set(CMAKE_REQUIRED_FLAGS "-I${OPENCL_INCLUDE_PATH}")
check_include_file_cxx(CL/cl.hpp      HAVE_OPENCL_H)
if(NOT HAVE_OPENCL_H)
  message(STATUS "WARNING: unable to find OpenCL header file.")
  message(STATUS "WARNING: OpenCL accelerated tools will NOT be built.")
endif(NOT HAVE_OPENCL_H)

## Basic library checks
include(CheckLibraryExists)

## C Math Library
set(MATH_LIBRARY_PATH "/usr/lib"
  CACHE STRING "Location of libm.so")
unset(HAVE_LIBM CACHE)
set(CMAKE_REQUIRED_FLAGS "-L${MATH_LIBRARY_PATH}")
check_library_exists(m sqrt "" HAVE_LIBM)
if(HAVE_LIBM)
  set(LIBS ${LIBS} m)
else(HAVE_LIBM)
  message(FATAL_ERROR "ERROR: unable to find C math library.")
endif(HAVE_LIBM)

## ARPACK
set(ARPACK_LIBRARY_PATH "/usr/lib"
  CACHE STRING "Location of libarpack.so")
unset(HAVE_LIBARPACK CACHE)
set(CMAKE_REQUIRED_FLAGS "-L${ARPACK_LIBRARY_PATH}")
check_library_exists(arpack dsaupd_ "" HAVE_LIBARPACK)
if(HAVE_LIBARPACK)
  set(LIBS ${LIBS} arpack)
else(HAVE_LIBARPACK)
  message(FATAL_ERROR "ERROR: unable to find ARPACK library.")
endif(HAVE_LIBARPACK)

## Oracle Berkeley DB
set(DB_CXX_LIBRARY_PATH "/usr/lib"
  CACHE STRING "Location of libdb_cxx.so")
set(CMAKE_REQUIRED_FLAGS "-L${DB_CXX_LIBRARY_PATH}")
unset(HAVE_LIBDB_CXX CACHE)
check_library_exists(db_cxx open "" HAVE_LIBDB_CXX)
if(HAVE_LIBDB_CXX)
  set(LIBS ${LIBS} db_cxx)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${DB_CXX_LIBRARY_PATH}")
else(HAVE_LIBDB_CXX)
  message(FATAL_ERROR "ERROR: unable to find Berkeley DB C++ library.")
endif(HAVE_LIBDB_CXX)

## OpenCL
set(OPENCL_LIBRARY_PATH "/usr/lib"
  CACHE STRING "Location of libOpenCL.so")
set(CMAKE_REQUIRED_FLAGS "-L${OPENCL_LIBRARY_PATH}")
unset(HAVE_LIBOPENCL CACHE)
if(HAVE_OPENCL_H)
check_library_exists(OpenCL clGetPlatformIDs "" HAVE_LIBOPENCL)
if(HAVE_LIBOPENCL)
  set(LIBS ${LIBS} OpenCL)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${OPENCL_LIBRARY_PATH}")
else(HAVE_LIBOPENCL)
  message(STATUS "WARNING: unable to find OpenCL library (libOpenCL.so).")
  message(STATUS "WARNING: OpenCL accelerated tools will NOT be built.")
endif(HAVE_LIBOPENCL)
endif(HAVE_OPENCL_H)

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

# GSL
pkg_check_modules(PC_GSL REQUIRED gsl)
set(GSL_DEFINITIONS ${PC_GSL_CFLAGS_OTHER})
find_path(GSL_INCLUDE_DIR gsl/gsl_vector.h
  HINTS ${PC_GSL_INCLUDEDIR} ${PC_GSL_INCLUDE_DIRS}
  PATH_SUFFIXES gromacs)
find_library(GSL_LIBRARY NAMES gsl libgsl
  HINTS ${PC_GSL_LIBDIR} ${PC_GSL_LIBRARY_DIRS})
find_library(GSL_CBLAS_LIBRARY NAMES gslcblas libgslcblas
  HINTS ${PC_GSL_LIBDIR} ${PC_GSL_LIBRARY_DIRS})
set(GSL_INCLUDE_DIRS ${GSL_INCLUDE_DIR})
set(GSL_LIBRARIES ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})

find_package_handle_standard_args(gsl DEFAULT_MSG
  GSL_LIBRARY GSL_INCLUDE_DIR)
mark_as_advanced(GSL_INCLUDE_DIR GSL_LIBRARY)

include_directories(${GSL_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GSL_LIBRARIES})

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

# Build tools and scripts...
set(MDSCTK_TOOLS
  auto_decomp_sparse
  auto_decomp_sparse_nystrom
  bb_xtc_to_phipsi
  check_xtc
  decomp_sparse
  decomp_sparse_nystrom
  knn_data
  knn_rms
  make_sysparse
  make_gesparse
  phipsi_to_sincos
  rms_test
  split_xtc)

foreach(TOOL ${MDSCTK_TOOLS})
  add_executable(${TOOL} ${TOOL}.cpp)
  target_link_libraries(${TOOL} ${LIBS})
endforeach(TOOL ${MDSCTK_TOOLS})

## OpenCL tools...
# set(MDSCTK_OCL_TOOLS
#   knn_data_ocl)

# foreach(TOOL ${MDSCTK_OCL_TOOLS})
#   add_executable(${TOOL} ${TOOL}.cpp ocl.o)
#   target_link_libraries(${TOOL} ${LIBS})
# endforeach(TOOL ${MDSCTK_OCL_TOOLS})
