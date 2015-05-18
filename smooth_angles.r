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
## Copyright (c) 2012-2015, Joshua L. Phillips.
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
    program.name <- "smooth_angles.r"
    source(paste(Sys.getenv("MDSCTK_HOME"),"/mdsctk.r",sep=""))
}

cat("   Perfoems amplitude- or frequency-based DFT smoothing for the\n")
cat("   set of angles in phipsi.dat. The size of the vectors is\n")
cat("   required. Amplitude-based smoothing is performed by default.\n")
cat("\n")
cat("   Use -h or --help to see the complete list of options.\n")
cat("\n")

parser$add_argument("-v","--vectorsize",type="integer",
                    help="Length of angle vectors",metavar="integer")
parser$add_argument("-p","--percent",default=10,type="double",
                    help="Smoothing percentage [default %(default)d]",metavar="double")
parser$add_argument("-f","--frequency", action="store_true", default=FALSE,
                    help="Use frequency-based smoothing (amplitude-based is the default)")
parser$add_argument("-a","--angles",default="phipsi.dat",
                    help="Angles file [default %(default)s]")
parser$add_argument("-o","--output",default="smoothedangles.dat",
                    help="(Output) Estimates file [default %(default)s]")

myargs <- parser$parse_args()

if (is.null(myargs$vectorsize)) {
    cat("ERROR: --vectorsize not supplied.\n")
    cat("\n")
    q()
}

myfrac <- myargs$percent / 100

cat("Running with the following options:\n")
cat(paste("vectorsize = ",myargs$vectorsize,"\n"))
if (myargs$frequency) {
    cat(paste("smoothing  = ","frequency-based","\n"))
} else {
    cat(paste("smoothing  = ","amplitude-based","\n"))
}
cat(paste("percent    = ",myargs$percent,"\n"))
cat(paste("angles     = ",myargs$angles,"\n"))
cat(paste("output     = ",myargs$output,"\n"))
cat("\n")

data <- apply(matrix(complex(modulus=1,argument=read.binary(myargs$angles,file.info(myargs$angles)$size/8)),nrow=myargs$vectorsize),1,fft)

if (myargs$frequency) {
    ## Frequency-based
    n <- round(nrow(data) * (1-myfrac) / 2)
    filter <- rep(0,nrow(data))
    filter[1] <- 1
    if (n > 0) {
        filter[seq(1,n)] <- 1
        filter[seq(nrow(data),nrow(data)-n+1)] <- 1
    }
    data <- Arg(apply(sweep(data,1,filter,FUN="*"),2,fft,inverse=TRUE)/nrow(data))
} else {
    ## Amplitude-based
    data.cutoff <- max(abs(data[-1,])) * myfrac
    data.dc <- data[1,]
    data[abs(data)<data.cutoff] <- 0
    data[1,] <- data.dc
    data <- Arg(apply(data,2,fft,inverse=TRUE)/nrow(data))
}

data[data < 0] <- data[data < 0] + 2*pi
myfd <- file(myargs$output,"wb")
writeBin(as.double(t(data)),myfd,size=8)
close(myfd)
