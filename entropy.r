#!/usr/bin/env Rscript
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
    program.name <- "entropy.r"
    source(paste(Sys.getenv("MDSCTK_HOME"),"/mdsctk.r",sep=""))
}

cat("   Computes the local entropy of the given sparse\n")
cat("   matrix with indices from indices.dat and the densities\n")
cat("   in density.dat. The number of nearest neighbors,\n")
cat("   k, is required.\n")
cat("\n")
cat("   Use -h or --help to see the complete list of options.\n")
cat("\n")

parser$add_argument("-k","--knn",type="integer",
                    help="Number of k-nearest neighbors",metavar="integer")
parser$add_argument("-i","--indices",default="indices.dat",
                    help="K-nn indices file [default %(default)s]")
parser$add_argument("-d","--densities",default="density.dat",
                    help="Density file [default %(default)s]")
parser$add_argument("-o","--output",default="entropy.dat",
                    help="(Output) Entropy file [default %(default)s]")

myargs <- parser$parse_args()

if (is.null(myargs$knn)) {
    cat("ERROR: --knn not supplied.\n")
    cat("\n")
    q()
}

cat("Running with the following options:\n")
cat(paste("knn =       ",myargs$knn,"\n"))
cat(paste("indices   = ",myargs$indices,"\n"))
cat(paste("densities = ",myargs$densities,"\n"))
cat(paste("output    = ",myargs$output,"\n"))
cat("\n")

myindexfile <- myargs$indices
mydensfile <- myargs$densities
myk <- myargs$knn

dens <- scan(mydensfile,quiet=TRUE)
write(entropy.pointwise(matrix(read.binary.int(myindexfile,length(dens)*myk)+1,ncol=myk,byrow=TRUE),dens),file=myargs$output,ncolumns=1)
