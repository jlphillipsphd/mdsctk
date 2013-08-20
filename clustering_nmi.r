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

if (Sys.getenv("MDSCTK_HOME")=="") {
    cat("\n")
    cat("Please set the MDSCTK_HOME environment variable\n")
    cat("before running this script.\n")
    cat("\n")
    q()
} else {
    source(paste(Sys.getenv("MDSCTK_HOME"),"/config.r",sep=""))
}

normalmutualinf <- function(data) {
  s <- 0
  e <- 0
  pc <- colSums(data)
  pr <- rowSums(data)
  for (i in 1:dim(data)[1])
    for (j in 1:dim(data)[2])
      if (data[i,j] > 0) {
        e <- e + (data[i,j] * log2(data[i,j]))
        s <- s + (data[i,j] * log2(data[i,j] / (pr[i] * pc[j])))
      }
  return(s / -e)
}

## myargs <- commandArgs(TRUE)

## if (length(myargs) != 1) {
##   print("")
##   print("Usage: clustering_nmi.r [cluster_histogram]")
##   print("")
##   q()
## }

## hist <- as.matrix(read.table(myargs[1],header=FALSE))
hist <- as.matrix(read.table("histogram.dat",header=FALSE))

myfd <- pipe("cat","wb")
write(normalmutualinf(hist),myfd,ncolumns=1)
close(myfd)
