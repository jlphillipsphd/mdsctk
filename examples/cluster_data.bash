#!/bin/bash
##
## 
##                This source code is part of
## 
##                        M D S C T K
## 
##       Molecular Dynamics Spectral Clustering ToolKit
## 
##                        VERSION 1.1.1
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

NTHREADS=2    ## Number of threads to use
KNN=20        ## Number of nearest neighbors to keep
NCLUSTERS=2   ## Number of clusters to extract
SCALING=10    ## (must be <= KNN) for calculating scaling factors...
DIM=2         ## Data dimensionality

echo "Computing distances between all point pairs..."
../knn_data ${NTHREADS} ${KNN} ${DIM} rings.pts rings.pts

echo "Creating CSC format sparse matrix..."
../make_sysparse ${KNN}

echo "Performing autoscaled spectral decomposition..."
../auto_decomp_sparse ${NCLUSTERS} ${SCALING}

echo "Clustering eigenvectors..."
../kmeans.r

## Make a plot of the results...
Rscript \
    -e 'data <- matrix(readBin("rings.pts",double(0),n=400*2),ncol=2,byrow=TRUE)' \
    -e 'clusters <- scan("clusters.dat")' \
    -e 'postscript("clusters.eps",width=5,height=5,onefile=FALSE,horizontal=FALSE)' \
    -e 'plot(data,col=c("red","blue")[clusters],pch=20)' \
    -e 'dev.off()'

echo "See clusters.eps for results (eg. evince clusters.eps)..."