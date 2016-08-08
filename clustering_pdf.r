#!/usr/bin/env Rscript
##
## 
##                This source code is part of
## 
##                        M D S C T K
## 
##       Molecular Dynamics Spectral Clustering ToolKit
## 
##                        VERSION 1.2.3
## Written by Joshua L. Phillips.
## Copyright (c) 2012-2016, Joshua L. Phillips.
## Check out http://www.cs.mtsu.edu/~jphillips/software.html for more
## information.
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
## For more info, check our website at
## http://www.cs.mtsu.edu/~jphillips/software.html
## 
##

if (Sys.getenv("MDSCTK_HOME")=="") {
    cat("\n")
    cat("Please set the MDSCTK_HOME environment variable\n")
    cat("before running this script.\n")
    cat("\n")
    q()
} else {
    program.name <- "clustering_pdf.r"
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
parser$add_argument("-o","--output",default="pdf.dat",
                    help="(Output) Joint probability distribution file [default %(default)s]")

myargs <- parser$parse_args()

cat("Running with the following options:\n")
cat(paste("assignment = ",myargs$assignment,"\n"))
cat(paste("clusters =   ",myargs$clusters,"\n"))
cat(paste("output =     ",myargs$output,"\n"))
cat("\n")

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
