#!/bin/bash
##
## 
##                This source code is part of
## 
##                        M D S C T K
## 
##       Molecular Dynamics Spectral Clustering ToolKit
## 
##                        VERSION 1.2.2
## Written by Joshua L. Phillips.
## Copyright (c) 2012-2014, Joshua L. Phillips.
## check out http://github.com/jlphillipsphd/mdsctk/ for more information.
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
## 
## If you want to redistribute modifications, please consider that
## derived work must not be called official MDSCTK. Details are found
## in the README & LICENSE files - if they are missing, get the
## official version at <website>.
## 
## To help us fund MDSCTK development, we humbly ask that you cite
## the papers on the package - you can find them in the top README file.
## 
## For more info, check our website at http://github.com/jlphillipsphd/mdsctk/
## 
##

if [[ -z "${MDSCTK_HOME}" ]]; then
    echo
    echo "Please set the MDSCTK_HOME environment variable"
    echo "before running this example..."
    echo
    exit 1
fi  

NTHREADS=8    ## Number of threads to use
KNN=1024      ## Number of nearest neighbors to keep
LENGTH=16     ## Length of the protein (number of AAs)
NANGLES=$((2*${LENGTH}-2))
DIM=$((${NANGLES}*2))  ## Data dimensionality
              ## (2*LENGTH-2)*2 (sin-cos of phi-psi angles)
WINDOW=1      ## Smoothing window for MLD
SMOOTH=50     ## Percentage of frequency smoothing

echo "## MLD Example from Phillips et al. 2014       ##"
echo "## This example will take some time to compute ##"
echo "## since the size of the data vectors is small ##"
echo "## which makes scaling rather difficult.       ##"
echo
echo "Figure 11 - GB1 dimensionality estimation results"
echo

echo "Unsmoothed results..."
echo "Converting protein frames into phi-psi angle space..."
bb_xtc_to_phipsi \
    -x traj/proteins/gb1.09.bb.xtc \
    -o gb1-unsmoothed-phipsi.dat
    
echo "Converting phi-psi angles to Euclidean sin-cos space.."
angles_to_sincos \
    -i gb1-unsmoothed-phipsi.dat \
    -o gb1-unsmoothed-sincos.dat

echo "Computing distances between all point pairs..."
#${MDSCTK_HOME}/knn_data_ocl -k ${KNN} -v ${DIM} \
${MDSCTK_HOME}/knn_data -t ${NTHREADS} -k ${KNN} -v ${DIM} \
    -r gb1-unsmoothed-sincos.dat \
    -d gb1-unsmoothed-distances.dat \
    -i gb1-unsmoothed-indices.dat

echo "Computing MLD for default k-neighbors..."
${MDSCTK_HOME}/mld_estimate.r -k ${KNN} -w ${WINDOW} \
    -d gb1-unsmoothed-distances.dat \
    -o gb1-unsmoothed-estimates.dat

echo "Frequency smoothing..."
echo "Smoothing Level - Frequency cutoff percent - %{SMOOTH}\%"

echo "Smoothing phi-psi angle space..."
smooth_angles.r -p ${SMOOTH} -f -v ${NANGLES} \
    -a gb1-unsmoothed-phipsi.dat \
    -o gb1-smoothed-phipsi.dat

echo "Converting phi-psi angles to Euclidean sin-cos space.."
angles_to_sincos \
    -i gb1-smoothed-phipsi.dat \
    -o gb1-smoothed-sincos.dat

echo "Computing distances between all point pairs..."
#${MDSCTK_HOME}/knn_data_ocl -k ${KNN} -v ${DIM} \
${MDSCTK_HOME}/knn_data -t ${NTHREADS} -k ${KNN} -v ${DIM} \
    -r gb1-smoothed-sincos.dat \
    -d gb1-smoothed-distances.dat \
    -i gb1-smoothed-indices.dat

echo "Computing MLD for default k-neighbors..."
${MDSCTK_HOME}/mld_estimate.r -k ${KNN} -w ${WINDOW} \
    -d gb1-smoothed-distances.dat \
    -o gb1-smoothed-estimates.dat

echo "Computing RMSD from main structure (folded-NFP or average-IDP)..."
${MDSCTK_HOME}/rms_test \
    -p traj/proteins/gb1.ca.gro \
    -x traj/proteins/gb1.09.ca.xtc \
    -o gb1-rmsd.dat

