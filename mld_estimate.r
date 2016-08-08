#!/usr/bin/env Rscript
##
## 
##                This source code is part of
## 
##                        M D S C T K
## 
##       Molecular Dynamics Spectral Clustering ToolKit
## 
##                        VERSION 1.2.5
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
## official version at github.com/jlphillipsphd/mdsctk/.
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
    program.name <- "mld_estimate.r"
    source(paste(Sys.getenv("MDSCTK_HOME"),"/mdsctk.r",sep=""))
}

cat("   Computes the Maximum Likelihood Dimension estimates for the\n")
cat("   sparse distances in distances.dat. The number of k nearest\n")
cat("   neighbors is required [knn]. The global MLD will be computed\n")
cat("   for each of the values in the [set] file. However, a default\n")
cat("   spread of 2,3,4,6,8,16,32,64,128,256,512,1024 is performed\n")
cat("   if the [set] file is not provided. Pointwise (local) estimates\n")
cat("   can be obtained by setting [window] to 1, and these estimates\n")
cat("   can be smoothed by averaging across several points by setting\n")
cat("   [window] > 1.\n")
cat("\n")
cat("   Use -h or --help to see the complete list of options.\n")
cat("\n")

parser$add_argument("-k","--knn",type="integer",
                    help="Number of k-nearest neighbors",metavar="integer")
parser$add_argument("-w","--window",default=0,type="integer",
                    help="Smoothing window [default %(default)d]",metavar="integer")
parser$add_argument("-s","--set",
                    help="Nearest neighbors file")
parser$add_argument("-d","--distances",default="distances.dat",
                    help="K-nn distances file [default %(default)s]")
parser$add_argument("-o","--output",default="mldestimates.dat",
                    help="(Output) Estimates file [default %(default)s]")

myargs <- parser$parse_args()

if (is.null(myargs$knn)) {
    cat("ERROR: --knn not supplied.\n")
    cat("\n")
    q()
}

cat("Running with the following options:\n")
cat(paste("knn       = ",myargs$knn,"\n"))
k <- c(2,3,4,6,8,16,32,64,128,256,512,1024)
if (!is.null(myargs$set))
    k <- scan(myargs$set,quiet=TRUE)
k <- k[k<=myargs$knn]
cat(paste("set       = ",paste(sprintf("%d",k),collapse=","),"\n"))
cat(paste("window    = ",myargs$window,"\n"))
cat(paste("distances = ",myargs$distances,"\n"))
cat(paste("output    = ",myargs$output,"\n"))
cat("\n")

myk <- myargs$knn

data <- matrix(read.binary(myargs$distances,file.info(myargs$distances)$size/8),nrow=myk)
result <- de.ml(data,k,window.size=myargs$window)
if (myargs$window > 0) {
    write(result,file=myargs$output,ncolumns=1)
} else {
    write(result,file=myargs$output,ncolumns=2)
}
