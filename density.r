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

myargs <- commandArgs(TRUE)

if (length(myargs) != 3 & length(myargs) != 4) {
  cat("\n")
  cat("Usage: density.r [distances.dat] [n] [k] <sigma>\n")
  cat("   Computes the kernel density estimate for the\n")
  cat("   provided sparse distance matrix.\n")
  cat("\n")
  q()
}

myfile <- myargs[1]
myn <- as.integer(myargs[2])
myk <- as.integer(myargs[3])

dens.pointwise <- function(data,sigma) {
  p <- colSums(exp(-data^2 / (2 * sigma^2)))
  return (p / nrow(data))

}

read.binary <- function(filename,n) {
  fd <- file(filename,open="rb")
  data <- readBin(fd,"double",n=n,size=8)
  close(fd)
  return (data)
}

data <- matrix(read.binary(myfile,myn*myk),ncol=myn)
if (length(myargs)==4) {
    sigma <- as.double(myargs[4])
}
else {
    sigma <- sd(as.double(data))
}
write(dens.pointwise(data,sigma),file=stdout(),ncolumns=1)
