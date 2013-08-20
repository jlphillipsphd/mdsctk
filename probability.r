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

if (length(myargs) != 1 & length(myargs) != 2) {
  cat("\n")
  cat(paste("   MDSCTK ",MDSCTK_VERSION_MAJOR,".",MDSCTK_VERSION_MINOR,"\n",sep=""))
  cat("   Copyright (C) 2013 Joshua L. Phillips\n")
  cat("   MDSCTK comes with ABSOLUTELY NO WARRANTY; see LICENSE for details.\n")
  cat("   This is free software, and you are welcome to redistribute it\n")
  cat("   under certain conditions; see README.md for details.\n")
  cat("\n")
  cat("Usage: density.r [k] <sigma>\n")
  cat("   Computes the kernel density estimate for the sparse distances\n")
  cat("   in distances.dat, and convert the results to estimated probailities.\n")
  cat("   The number of k nearest neighbors is required but the bandwidth of\n")
  cat("   the kernel, sigma, can be supplied or guesstimated based the data.\n")
  cat("\n")
  q()
}

myk <- as.integer(myargs[1])

prob.pointwise <- function(data,sigma) {
  p <- colSums(exp(-data^2 / (2 * sigma^2)))
  return (p / sum(p))

}

read.binary <- function(filename,n) {
  fd <- file(filename,open="rb")
  data <- readBin(fd,"double",n=n,size=8)
  close(fd)
  return (data)
}

data <- matrix(read.binary("distances.dat",file.info("distances.dat")$size/8),nrow=myk)
if (length(myargs)==2) {
    sigma <- as.double(myargs[2])
} else {
    sigma <- sd(as.double(data))
}
write(prob.pointwise(data,sigma),file="probability.dat",ncolumns=1)
