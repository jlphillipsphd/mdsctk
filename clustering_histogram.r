#!/usr/bin/env Rscript
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
## official version at github.com/douradopalmares/mdsctk/.
## 
## To help us fund MDSCTK development, we humbly ask that you cite
## the papers on the package - you can find them in the top README file.
## 
## For more info, check our website at http://github.com/douradopalmares/mdsctk/
## 
##

## Might utlize command-line arguments in future versions...
## myargs <- commandArgs(TRUE)

## if (length(myargs) != 2) {
##   print("")
##   print("Usage: clustering_histogram.r [traj_assignment] [cluster_assignment]")
##   print("")
##   q()
## }

## Utilities
cluster.sort <- function(data,f=median) {
  s <- seq(range(data)[1], range(data)[2])
  x <- rep(0,length(s))
  l <- list()
  for (i in 1:length(s)) {
    l[[i]] <- which(data == s[i])
    x[i] <- f(l[[i]])
  }
  ix <- sort(x, index.return=TRUE)$ix
  for (i in 1:length(s)) data[l[[ix[i]]]] <- i
  return (data)
}

## clus.assign <- cluster.sort(scan(myargs[2],quiet=TRUE))
## myassignment <- scan(myargs[1],quiet=TRUE)

## Note that the raw cluster numbers are changed, or "resorted"
## in order to make the histogram plots interpretable. These can
## be kept with their original labels by simply not using the
## cluster.sort() function, which will not change the NMI of the
## result, but will make the plots more difficult to interpret.
clus.assign <- cluster.sort(scan("clusters.dat",quiet=TRUE))
myassignment <- scan("assignment.dat",quiet=TRUE)
nclusters <- max(clus.assign)
ntraj <- max(myassignment)
nframes <- length(clus.assign)

myhist <- matrix(0,nclusters,ntraj)
for (x in seq(1,length(myassignment))) {
  myhist[clus.assign[x],myassignment[x]] <- myhist[clus.assign[x],myassignment[x]] + 1 
}
myhist <- sweep(myhist,2,colSums(myhist),FUN="/")
myhist <- myhist / sum(myhist)

myfd <- pipe("cat","wb")
myfd <- file("histogram.dat","w")
# Note that the matrix is transposed when written due
# to column-order processing (R and Fortran).
write(myhist,myfd,ncolumns=nclusters)
close(myfd)
