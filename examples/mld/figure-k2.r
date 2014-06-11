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
## official version at <website>.
## 
## To help us fund MDSCTK development, we humbly ask that you cite
## the papers on the package - you can find them in the top README file.
## 
## For more info, check our website at http://github.com/douradopalmares/mdsctk/
## 
##

myargs <- commandArgs(TRUE)
unsmoothed.file <- myargs[1]
smoothed.file <- myargs[2]
figure.file <- myargs[3]

unsmoothed <- read.table(unsmoothed.file,header=FALSE)
smoothed <- read.table(smoothed.file,header=FALSE)

## Clean up data...
k <- unsmoothed[,1]

unsmoothed <- unsmoothed[,seq(2,ncol(unsmoothed),2)]
smoothed <- smoothed[,seq(2,ncol(smoothed),2)]

row.names(unsmoothed) <- k
row.names(smoothed) <- k
names(unsmoothed) <- c("0.0","0.01","0.1","1.0","3.0","10.0")
names(smoothed) <- c("0.0","0.01","0.05","0.1","0.5")

## Make plot
postscript(figure.file,width=3.35,height=5,onefile=FALSE,horizontal=FALSE)
layout(matrix(c(1,2)),heights=c(1,1.15))

## First plot
par(mar=c(0,2,0.1,0)+0.1,mgp=c(1.1,0.5,0),cex.axis=0.5,cex.lab=0.75)
matplot(as.double(row.names(unsmoothed)),unsmoothed,type='o',
        lty=seq(1,ncol(unsmoothed)),lwd=1,pch=seq(19,24),cex=0.5,
        xlab="k",
        ylab="Estimated Dimensionality",ylim=range(unsmoothed),
        xaxt="n")
legend("topright",names(unsmoothed),lty=seq(1,ncol(unsmoothed)),lwd=1,cex=0.5,col=seq(1,ncol(unsmoothed)),bg="white",pch=seq(19,24),title=expression(sigma[theta[unfolded]]))
mtext("A",outer=TRUE,adj=0,line=-1)

## Second plot
par(mar=c(2,2,0,0)+0.1,mgp=c(1.1,0.5,0),cex.axis=0.5,cex.lab=0.75)
matplot(as.double(row.names(smoothed)),smoothed,type='o',
        lty=seq(1,ncol(smoothed)),lwd=1,pch=seq(19,24),cex=0.5,
        xlab="k",
        ylab="Estimated Dimensionality",ylim=range(smoothed))
legend("topright",names(smoothed),lty=seq(1,ncol(smoothed)),lwd=1,cex=0.5,col=seq(1,ncol(smoothed)),bg="white",pch=seq(19,24),title="Fraction")
mtext("B",outer=TRUE,adj=0,line=-12.6)

dev.off()
