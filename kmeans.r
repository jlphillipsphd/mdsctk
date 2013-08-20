#!/usr/bin/env Rscript
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
## official version at github.com/douradopalmares/mdsctk/.
## 
## To help us fund MDSCTK development, we humbly ask that you cite
## the papers on the package - you can find them in the top README file.
## 
## For more info, check our website at http://github.com/douradopalmares/mdsctk/
## 
##

if (Sys.getenv("MDSCTK_HOME")=="") {
    cat("\n")
    cat("Please set the MDSCTK_HOME environment variable\n")
    cat("before running this script.\n")
    cat("\n")
    q()
} else {
    source(paste(Sys.getenv("MDSCTK_HOME"),"/config.r",sep=""))
}

myargs <- commandArgs(TRUE)

if (length(myargs) != 1) {
  cat("\n")
  cat(paste("   MDSCTK ",MDSCTK_VERSION_MAJOR,".",MDSCTK_VERSION_MINOR,"\n",sep=""))
  cat("   Copyright (C) 2013 Joshua L. Phillips\n")
  cat("   MDSCTK comes with ABSOLUTELY NO WARRANTY; see LICENSE for details.\n")
  cat("   This is free software, and you are welcome to redistribute it\n")
  cat("   under certain conditions; see README.md for details.\n")
  cat("\n")
  cat("Usage: kmeans.r [k]\n")
  cat("   Performs standard k-means clustering on the provided\n")
  cat("   eigenvectors from (eigenvalues.dat and eigenvectors.dat).\n")
  cat("   The number of clusters requested (k), can be 2>=k<=nev\n")
  cat("   where nev is the number of eigenvectors. The results are\n")
  cat("   written to clusters.dat, and a breakdown of assignments\n")
  cat("   by cluster is written to clusters.ndx.\n")
  cat("\n")
  q()
}

valin <- "eigenvalues.dat"
vecin <- "eigenvectors.dat"
nclusters <- as.integer(myargs[1])

e.values <- as.double(scan(valin,quiet=TRUE))

if (nclusters > length(e.values) | nclusters < 2) {
  cat("\n")
  cat("   -- ERROR-- \n")
  cat(sprintf("   Number of eigenvalues: %d\n",length(e.values)))
  cat(sprintf("   Number of clusters requested: %d\n",nclusters))
  cat("   The number of clusters requested must be >=2 and\n")
  cat("   <= number of eigenvalues.\n")
  cat("\n")
  q()    
}

e.vectors <- matrix(scan(vecin,quiet=TRUE),ncol=length(e.values))

set.seed(0) # Change for different results...
clusters <- kmeans(e.vectors[,seq(1,nclusters)],
                   nclusters,
                   iter.max=30,
                   nstart=10)$cluster

myout <- file("clusters.dat","w")
write(clusters,myout,ncolumns=1)
close(myout)

myout <- file("clusters.ndx","w")
for (n in seq(1,ncol(e.vectors))) {
  writeLines(sprintf("[cluster_%d]",n),con=myout)
  write(which(clusters==n),myout,ncolumns=20)
  writeLines("",con=myout)
}
close(myout)

q(save="no")
