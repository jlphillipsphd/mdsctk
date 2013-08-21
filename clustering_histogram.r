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
    program.name <- "clustering_histogram.r"
    source(paste(Sys.getenv("MDSCTK_HOME"),"/mdsctk.r",sep=""))
}

cat("   Reads in a set of assignments (classes,groups,etc.) for each\n")
cat("   data point or structure from an assignment file, and the cluster\n")
cat("   assignment data from a cluster assignment file in order to calculate\n")
cat("   a joint cluster-assignment probability distribution which\n")
cat("   is written to the provided output file.\n")
cat("   Note that cluster assignments are relabeled to better\n")
cat("   visualize the mutual information between assignments\n")
cat("   and clusters.\n")
cat("\n")
cat("   Use -h or --help to see the complete list of options.\n")
cat("\n")

parser$add_argument("-a","--assignment",default="assignment.dat",
                    help="Data assignment file [default %(default)s]")
parser$add_argument("-c","--clusters",default="clusters.dat",
                    help="Cluster assignment file [default %(default)s]")
parser$add_argument("-o","--output",default="histogram.dat",
                    help="(Output) Joint probability distribution file [default %(default)s]")

myargs <- parser$parse_args()

cat("Running with the following options:\n")
cat(paste("assignment = ",myargs$assignment,"\n"))
cat(paste("clusters =   ",myargs$clusters,"\n"))
cat(paste("output =     ",myargs$output,"\n"))
cat("\n")

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

## Note that the raw cluster numbers are "changed", or relabeled
## in order to make the histogram plots interpretable. These can
## be kept with their original labels by simply not using the
## cluster.sort() function, which will not change the NMI of the
## result, but will make the plots more difficult to interpret.
clus.assign <- cluster.sort(scan(myargs$clusters,quiet=TRUE))
myassignment <- scan(myargs$assignment,quiet=TRUE)
nclusters <- max(clus.assign)
ntraj <- max(myassignment)
nframes <- length(clus.assign)

myhist <- matrix(0,nclusters,ntraj)
for (x in seq(1,length(myassignment))) {
  myhist[clus.assign[x],myassignment[x]] <- myhist[clus.assign[x],myassignment[x]] + 1 
}
myhist <- sweep(myhist,2,colSums(myhist),FUN="/")
myhist <- myhist / sum(myhist)

myfd <- file(myargs$output,"w")
# Note that the matrix is transposed when written due
# to column-order processing (R and Fortran).
write(myhist,myfd,ncolumns=nclusters)
close(myfd)
