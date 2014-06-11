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
## check out http://github.com/douradopalmares/mdsctk/ for more information.
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
## For more info, check our website at http://github.com/douradopalmares/mdsctk/
## 
##

if [[ -z "${MDSCTK_HOME}" ]]; then
    echo
    echo "Please set the MDSCTK_HOME environment variable"
    echo "before running this example..."
    echo
    exit 1
fi  

NTHREADS=2    ## Number of threads to use
KNN=1024      ## Number of nearest neighbors to keep
POLYMER=20    ## Length (number of points) in the
              ## polymer (one more than number of links)
NANGLES=$((2*${POLYMER}-5))
DIM=$((${NANGLES}*2))  ## Data dimensionality
              ## (2*POLYMER-5)*2 (sin-cos of theta-phi angles)

echo "## MLD Example from Phillips et al. 2014       ##"
echo "## This example will take some time to compute ##"
echo "## since the size of the data vectors is small ##"
echo "## which makes scaling rather difficult.       ##"
echo
echo "Figure 4 - Half-folded helix model results"
echo

for NOISE in 00.00 00.01 00.10 01.00 03.00 10.00; do
    echo "Unsmoothed results..."
    echo "Noise Level - sigma_{theta,phi} = ${NOISE}"
    echo "Converting polymer frames into theta-phi angle space..."
    ca_xtc_to_thetaphi \
	-x traj/polymers/half_folded_helix/${POLYMER}/${NOISE}/chain.5000.xtc \
	-o halffolded-unsmoothed-${NOISE}-thetaphi.dat
    
    echo "Converting theta-phi angles to Euclidean sin-cos space.."
    angles_to_sincos \
	-i halffolded-unsmoothed-${NOISE}-thetaphi.dat \
	-o halffolded-unsmoothed-${NOISE}-sincos.dat

    echo "Computing distances between all point pairs..."
    #${MDSCTK_HOME}/knn_data_ocl -k ${KNN} -v ${DIM} \
    ${MDSCTK_HOME}/knn_data -t ${NTHREADS} -k ${KNN} -v ${DIM} \
	-r halffolded-unsmoothed-${NOISE}-sincos.dat \
	-d halffolded-unsmoothed-${NOISE}-distances.dat \
	-i halffolded-unsmoothed-${NOISE}-indices.dat

    echo "Computing MLD for default k-neighbors..."
    ${MDSCTK_HOME}/mld_estimate.r -k ${KNN} \
	-d halffolded-unsmoothed-${NOISE}-distances.dat \
	-o halffolded-unsmoothed-${NOISE}-estimates.dat
done

## Combine the results
paste halffolded-unsmoothed-??.??-estimates.dat > halffolded-unsmoothed.dat

NOISE=00.10
for SMOOTH in 00 01 05 10 50; do
    echo "Frequency smoothing..."
    echo "Smoothing Level - Frequency cutoff percent - ${SMOOTH}\% (sigma_{theta,phi} = ${NOISE})"

    echo "Smoothing theta-phi angle space..."
    smooth_angles.r -p ${SMOOTH} -f -v ${NANGLES} \
	-a halffolded-unsmoothed-${NOISE}-thetaphi.dat \
	-o halffolded-smoothed-${SMOOTH}-thetaphi.dat
    
    echo "Converting theta-phi angles to Euclidean sin-cos space.."
    angles_to_sincos \
	-i halffolded-smoothed-${SMOOTH}-thetaphi.dat \
	-o halffolded-smoothed-${SMOOTH}-sincos.dat

    echo "Computing distances between all point pairs..."
    #${MDSCTK_HOME}/knn_data_ocl -k ${KNN} -v ${DIM} \
    ${MDSCTK_HOME}/knn_data -t ${NTHREADS} -k ${KNN} -v ${DIM} \
	-r halffolded-smoothed-${SMOOTH}-sincos.dat \
	-d halffolded-smoothed-${SMOOTH}-distances.dat \
	-i halffolded-smoothed-${SMOOTH}-indices.dat

    echo "Computing MLD for default k-neighbors..."
    ${MDSCTK_HOME}/mld_estimate.r -k ${KNN} \
	-d halffolded-smoothed-${SMOOTH}-distances.dat \
	-o halffolded-smoothed-${SMOOTH}-estimates.dat
done

## Combine the results
paste halffolded-smoothed-??-estimates.dat > halffolded-smoothed.dat

