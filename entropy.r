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
  cat("Usage: entropy.r [k]\n")
  cat("   Computes the local entropy of the given sparse\n")
  cat("   matrix with indices from indices.dat and the densities\n")
  cat("   in density.dat. The number of nearest neighbors,\n")
  cat("   k, is required.\n")
  cat("\n")
  q()
}

myindexfile <- "indices.dat"
mydensfile <- "density.dat"
myk <- as.integer(myargs[1])

entropy.pointwise <- function(data,dens) {
  return (apply(apply(data,
                      1,
                      FUN=function(vec,dens){dens[vec]},
                      dens),
                2,FUN=function(vec){-sum(log2(vec))/length(vec)}))
}

read.binary.int <- function(filename,n) {
  fd <- file(filename,open="rb")
  data <- readBin(fd,"integer",n=n,size=4)
  close(fd)
  return (data)
}

dens <- scan(mydensfile,quiet=TRUE)
write(entropy.pointwise(matrix(read.binary.int(myindexfile,length(dens)*myk)+1,ncol=myk,byrow=TRUE),dens),file="entropy.dat",ncolumns=1)
