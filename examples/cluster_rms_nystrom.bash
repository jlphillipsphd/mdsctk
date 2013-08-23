#!/bin/bash
##
## 
##                This source code is part of
## 
##                        M D S C T K
## 
##       Molecular Dynamics Spectral Clustering ToolKit
## 
##                        VERSION 1.1.2
## Written by Joshua L. Phillips.
## Copyright (c) 2013, Joshua L. Phillips.
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

TOP=trp-cage.pdb
XTC=trp-cage.xtc
OSXTC=trp-cage-outofsample.xtc

NTHREADS=2    # The number of threads to use
KNN=100       # The number of nearest-neighbors to keep for
              # the sparse representation
NCLUSTERS=10  # The number of clusters to extract
SCALING=12    # The number of nearest-neighbors to use for
              # computing scaling factors (should be < KNN)

echo "Computing RMS distances between all landmark pairs..."
${MDSCTK_HOME}/knn_rms -t ${NTHREADS} -k ${KNN} -p ${TOP} -r ${XTC}

echo "Creating CSC format symmetric sparse matrix..."
${MDSCTK_HOME}/make_sysparse -k ${KNN}

echo "Computing RMS distances between landmarks and out-of-sample structures..."
${MDSCTK_HOME}/knn_rms -t ${NTHREADS} -k ${KNN} -p ${TOP} -r ${XTC} -f ${OSXTC}

echo "Creating CSC format non-symmetric sparse matrix..."
${MDSCTK_HOME}/make_gesparse -k ${KNN}

echo "Performing autoscaled spectral decomposition..."
${MDSCTK_HOME}/auto_decomp_sparse_nystrom -n ${NCLUSTERS} -k ${SCALING}

# echo "Performing spectral decomposition..."
# ${MDSCTK_HOME}/decomp_sparse_nystrom -n ${NCLUSTERS} -q 5.0

echo "Clustering eigenvectors..."
${MDSCTK_HOME}/kmeans.r -k ${NCLUSTERS}

# Generate trajectory assignment file,
# 10 trajectories of 100 frames each
# and a single remaining out-of-sample
# trajectory with 10000 frames.
Rscript \
    -e 'myout<-pipe("cat","w")' \
    -e 'write(c(sort(rep(seq(1,10),100)),rep(11,10000)),myout,ncolumns=1)' \
    -e 'close(myout)' > assignment.dat

echo "Computing replicate-cluster assignment histogram..."
${MDSCTK_HOME}/clustering_histogram.r

echo "Plotting the histogram (fails if R package 'fields' is missing)..."
${MDSCTK_HOME}/plot_histogram.r

echo "See histogram.eps for results (eg. evince histogram.eps)..."

echo "Computing normalized mutual information..."
${MDSCTK_HOME}/clustering_nmi.r | tee nmi.dat
