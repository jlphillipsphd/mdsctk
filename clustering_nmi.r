#!/usr/bin/env Rscript
##
## 
##                This source code is part of
## 
##                        M D S C T K
## 
##       Molecular Dynamics Spectral Clustering ToolKit
## 
##                        VERSION 1.2.2
## Written by Joshua L. Phillips.
## Copyright (c) 2012-2014, Joshua L. Phillips.
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
    program.name <- "clustering_nmi.r"
    source(paste(Sys.getenv("MDSCTK_HOME"),"/mdsctk.r",sep=""))
}

cat("   Reads in a joint assignment-cluster probability distribution\n")
cat("   from and computes the normalized mutual information between\n")
cat("   assignments and clusters.\n")
cat("\n")
cat("   Use -h or --help to see the complete list of options.\n")
cat("\n")

parser$add_argument("-p","--pdf",default="pdf.dat",
                    help="Joint PDF file [default %(default)s]")

myargs <- parser$parse_args()

cat("Running with the following options:\n")
cat(paste("pdf = ",myargs$pdf,"\n"))
cat("\n")

hist <- as.matrix(read.table(myargs$pdf,header=FALSE))

cat(sprintf("NMI = %g\n",normalmutualinf(hist)))
cat("\n")
