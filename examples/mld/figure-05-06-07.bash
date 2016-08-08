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

NTHREADS=2    ## Number of threads to use
KNN=1024      ## Number of nearest neighbors to keep
POLYMER=20    ## Length (number of points) in the
              ## polymer (one more than number of links)
NANGLES=$((2*${POLYMER}-5))
DIM=$((${NANGLES}*2))  ## Data dimensionality
              ## (2*POLYMER-5)*2 (sin-cos of theta-phi angles)
WINDOW=500    ## Smoothing window size

echo "## MLD Example from Phillips et al. 2014       ##"
echo "## This example will take some time to compute ##"
echo "## since the size of the data vectors is small ##"
echo "## which makes scaling rather difficult.       ##"
echo
echo "Figure 5 - Correlated helix model results"
echo "Figure 6 - Correlated helix model results (pointwise-unsmoothed)"
echo "Figure 7 - Correlated helix model results (pointwise-smoothed)"
echo

for NOISE in 00.00 00.01 00.10 01.00 03.00 10.00; do
    echo "Unsmoothed results..."
    echo "Noise Level - sigma_{theta,phi} = ${NOISE}"
    echo "Converting polymer frames into theta-phi angle space..."
    ca_xtc_to_thetaphi \
    	-x traj/polymers/correlated_helix/${POLYMER}_5/${NOISE}/chain.5000.xtc \
    	-o correlated-unsmoothed-${NOISE}-thetaphi.dat
    
    echo "Converting theta-phi angles to Euclidean sin-cos space.."
    angles_to_sincos \
    	-i correlated-unsmoothed-${NOISE}-thetaphi.dat \
    	-o correlated-unsmoothed-${NOISE}-sincos.dat

    echo "Computing distances between all point pairs..."
    #${MDSCTK_HOME}/knn_data_ocl -k ${KNN} -v ${DIM} \
    ${MDSCTK_HOME}/knn_data -t ${NTHREADS} -k ${KNN} -v ${DIM} \
    	-r correlated-unsmoothed-${NOISE}-sincos.dat \
    	-d correlated-unsmoothed-${NOISE}-distances.dat \
    	-i correlated-unsmoothed-${NOISE}-indices.dat

    echo "Computing MLD for default k-neighbors..."
    ${MDSCTK_HOME}/mld_estimate.r -k ${KNN} \
    	-d correlated-unsmoothed-${NOISE}-distances.dat \
    	-o correlated-unsmoothed-${NOISE}-estimates.dat

    echo "Computing pointwise MLD for small-k..."
    ${MDSCTK_HOME}/mld_estimate.r -k ${KNN} -s small_k.txt -w ${WINDOW} \
	-d correlated-unsmoothed-${NOISE}-distances.dat \
	-o correlated-unsmoothed-${NOISE}-estimates-small.dat

    echo "Computing pointwise MLD for medium-k..."
    ${MDSCTK_HOME}/mld_estimate.r -k ${KNN} -s medium_k.txt -w ${WINDOW} \
	-d correlated-unsmoothed-${NOISE}-distances.dat \
	-o correlated-unsmoothed-${NOISE}-estimates-medium.dat

    echo "Computing pointwise MLD for large-k..."
    ${MDSCTK_HOME}/mld_estimate.r -k ${KNN} -s large_k.txt -w ${WINDOW} \
	-d correlated-unsmoothed-${NOISE}-distances.dat \
	-o correlated-unsmoothed-${NOISE}-estimates-large.dat
done

## Combine the results
paste correlated-unsmoothed-??.??-estimates.dat > correlated-unsmoothed.dat
paste correlated-unsmoothed-??.??-estimates-small.dat > correlated-unsmoothed-small.dat
paste correlated-unsmoothed-??.??-estimates-medium.dat > correlated-unsmoothed-medium.dat
paste correlated-unsmoothed-??.??-estimates-large.dat > correlated-unsmoothed-large.dat

NOISE=00.10
for SMOOTH in 00 01 05 10 50; do
    echo "Frequency smoothing..."
    echo "Smoothing Level - Frequency cutoff percent - ${SMOOTH}\% (sigma_{theta,phi} = ${NOISE})"

    echo "Smoothing theta-phi angle space..."
    smooth_angles.r -p ${SMOOTH} -f -v ${NANGLES} \
    	-a correlated-unsmoothed-${NOISE}-thetaphi.dat \
    	-o correlated-smoothed-${SMOOTH}-thetaphi.dat
    
    echo "Converting theta-phi angles to Euclidean sin-cos space.."
    angles_to_sincos \
    	-i correlated-smoothed-${SMOOTH}-thetaphi.dat \
    	-o correlated-smoothed-${SMOOTH}-sincos.dat

    echo "Computing distances between all point pairs..."
    #${MDSCTK_HOME}/knn_data_ocl -k ${KNN} -v ${DIM} \
    ${MDSCTK_HOME}/knn_data -t ${NTHREADS} -k ${KNN} -v ${DIM} \
    	-r correlated-smoothed-${SMOOTH}-sincos.dat \
    	-d correlated-smoothed-${SMOOTH}-distances.dat \
    	-i correlated-smoothed-${SMOOTH}-indices.dat

    echo "Computing MLD for default k-neighbors..."
    ${MDSCTK_HOME}/mld_estimate.r -k ${KNN} \
    	-d correlated-smoothed-${SMOOTH}-distances.dat \
    	-o correlated-smoothed-${SMOOTH}-estimates.dat

    echo "Computing pointwise MLD for small-k..."
    ${MDSCTK_HOME}/mld_estimate.r -k ${KNN} -s small_k.txt -w ${WINDOW} \
	-d correlated-smoothed-${SMOOTH}-distances.dat \
	-o correlated-smoothed-${SMOOTH}-estimates-small.dat

    echo "Computing pointwise MLD for medium-k..."
    ${MDSCTK_HOME}/mld_estimate.r -k ${KNN} -s medium_k.txt -w ${WINDOW} \
	-d correlated-smoothed-${SMOOTH}-distances.dat \
	-o correlated-smoothed-${SMOOTH}-estimates-medium.dat

    echo "Computing pointwise MLD for large-k..."
    ${MDSCTK_HOME}/mld_estimate.r -k ${KNN} -s large_k.txt -w ${WINDOW} \
	-d correlated-smoothed-${SMOOTH}-distances.dat \
	-o correlated-smoothed-${SMOOTH}-estimates-large.dat    
done

## Combine the results
paste correlated-smoothed-??-estimates.dat > correlated-smoothed.dat
paste correlated-smoothed-??-estimates-small.dat > correlated-smoothed-small.dat
paste correlated-smoothed-??-estimates-medium.dat > correlated-smoothed-medium.dat
paste correlated-smoothed-??-estimates-large.dat > correlated-smoothed-large.dat

