#!/bin/bash
##
## 
##                This source code is part of
## 
##                        M D S C T K
## 
##       Molecular Dynamics Spectral Clustering ToolKit
## 
##                        VERSION 1.2.0
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

NTHREADS=2    ## Number of threads to use
KNN=12        ## Number of nearest neighbors to keep
DIM=3         ## Data dimensionality

echo "Computing distances between all point pairs..."
${MDSCTK_HOME}/knn_data -t ${NTHREADS} -k ${KNN} -v ${DIM} -r swissroll.pts

echo "Creating CSC format sparse matrix..."
${MDSCTK_HOME}/make_sysparse -k ${KNN}

echo "Calculating approximate geodesic distances..."
${MDSCTK_HOME}/dijkstra -t ${NTHREADS}

echo "Performing spectral decomposition..."
${MDSCTK_HOME}/decomp_dense

## Make a plot of the resulting 2D embedding...
Rscript \
    -e 'data <- matrix(scan("eigenvectors.dat",n=1000*2),ncol=2)' \
    -e 'postscript("isomap.eps",width=5,height=5,onefile=FALSE,horizontal=FALSE)' \
    -e 'plot(data,col=rainbow(1000),pch=5)' \
    -e 'temp <- dev.off()'

echo "See isomap.eps for results (eg. evince isomap.eps)..."