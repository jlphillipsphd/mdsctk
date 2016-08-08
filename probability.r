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
    program.name <- "probability.r"
    source(paste(Sys.getenv("MDSCTK_HOME"),"/mdsctk.r",sep=""))
}

cat("   Computes the kernel density estimate for the sparse distances\n")
cat("   in distances.dat, and convert the results to estimated probailities.\n")
cat("   The number of k nearest neighbors is required but the bandwidth of\n")
cat("   the kernel, sigma, can be supplied or guesstimated based the data.\n")
cat("\n")
cat("   Use -h or --help to see the complete list of options.\n")
cat("\n")

parser$add_argument("-k","--knn",type="integer",
                    help="Number of k-nearest neighbors",metavar="integer")
parser$add_argument("-q","--sigma",type="double",
                    help="Kernel bandwidth",metavar="number")
parser$add_argument("-i","--indices",default="indices.dat",
                    help="K-nn indices file [default %(default)s]")
parser$add_argument("-d","--distances",default="distances.dat",
                    help="K-nn distances file [default %(default)s]")
parser$add_argument("-o","--output",default="probability.dat",
                    help="(Output) Probability file [default %(default)s]")

myargs <- parser$parse_args()

if (is.null(myargs$knn)) {
    cat("ERROR: --knn not supplied.\n")
    cat("\n")
    q()
}

cat("Running with the following options:\n")
cat(paste("knn =       ",myargs$knn,"\n"))
if (is.null(myargs$sigma)) {
    cat(paste("sigma =     ","estimated","\n"))
} else {
    cat(paste("sigma =     ",myargs$sigma,"\n"))
}
cat(paste("indices   = ",myargs$indices,"\n"))
cat(paste("distances = ",myargs$distances,"\n"))
cat(paste("output    = ",myargs$output,"\n"))
cat("\n")

myk <- myargs$knn

data <- matrix(read.binary(myargs$distances,file.info(myargs$distances)$size/8),nrow=myk)
if (is.null(myargs$sigma)) {
    sigma <- sd(as.double(data))
    cat(sprintf("Estimated sigma: %g\n\n",sigma))
} else {
    sigma <- myargs$sigma
}
write(prob.pointwise(data,sigma),file=myargs$output,ncolumns=1)
