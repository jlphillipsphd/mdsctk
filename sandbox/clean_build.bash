#!/bin/bash

## Scripts and text files
rm MDSCTK.bash
rm config.h
rm config.r
rm Makefile

## Executables
rm angles_to_sincos \
   auto_decomp_sparse \
   auto_decomp_sparse_nystrom \
   auto_heir_decomp_sparse \
   bb_xtc_to_phipsi \
   ca_xtc_to_thetaphi \
   check_xtc \
   contact_profile \
   decomp_dense \
   decomp_sparse \
   decomp_sparse_nystrom \
   dijkstra \
   flatten_xtc \
   knn_data \
   knn_data_sparse \
   knn_rms \
   make_sysparse \
   make_gesparse \
   rms_test \
   split_xtc

## CMake files
rm CMakeCache.txt
rm -rf CMakeFiles
rm cmake_install.cmake
