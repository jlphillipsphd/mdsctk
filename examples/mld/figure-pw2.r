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
small.file <- myargs[1]
medium.file <- myargs[2]
large.file <- myargs[3]
window <- as.integer(myargs[4])
figure.file <- myargs[5]

small <- read.table(small.file,header=FALSE)
medium <- read.table(medium.file,header=FALSE)
large <- read.table(large.file,header=FALSE)

## Clean up data...
names(small) <- c("0.0","0.01","0.05","0.1","0.5")
names(medium) <- c("0.0","0.01","0.05","0.1","0.5")
names(large) <- c("0.0","0.01","0.05","0.1","0.5")

row.names(small) <- seq(window,5000)
row.names(medium) <- seq(window,5000)
row.names(large) <- seq(window,5000)

mymin <- min(c(min(small),min(medium),min(large)))
mymax <- max(c(max(small),max(medium),max(large)))

## Make plot
postscript(figure.file,width=3.35,height=5.8,onefile=FALSE,horizontal=FALSE)
layout(matrix(c(1,2,3)),heights=c(1.5,1.5,1.7))

## First plot
par(mar=c(0,2.3,0,0)+0.1,mgp=c(1.1,0.5,0),cex.axis=0.75,cex.lab=1)
matplot(seq(500,5000),small,type='l',lwd=1,xlab="Time Step (t)",
        ylab="Estimated Dimensionality",ylim=c(mymin,mymax),lty=seq(1,ncol(small)),
        xaxt="n")
## legend("bottomright",names(small),lty=seq(1,ncol(small)),lwd=1,cex=0.5,
##        col=seq(1,ncol(small)),bg="white",title="Fraction")
mtext("A",outer=TRUE,adj=0,line=-1.2)

## Second plot
par(mar=c(0,2.3,0,0)+0.1,mgp=c(1.1,0.5,0),cex.axis=0.75,cex.lab=1)
matplot(seq(500,5000),medium,type='l',lwd=1,xlab="Time Step (t)",
        ylab="Estimated Dimensionality",ylim=c(mymin,mymax),lty=seq(1,ncol(medium)),
        xaxt="n")
## legend("bottomright",names(medium),lty=seq(1,ncol(medium)),lwd=1,
##        col=seq(1,ncol(medium)),bg="white",title="Fraction")
mtext("B",outer=TRUE,adj=0,line=-15.4)

## Third plot
par(mar=c(2,2.3,0,0)+0.1,mgp=c(1.1,0.5,0),cex.axis=0.75,cex.lab=1)
matplot(seq(500,5000),large,type='l',lwd=1,xlab="Time Step (t)",
        ylab="Estimated Dimensionality",ylim=c(mymin,mymax),lty=seq(1,ncol(large)))
legend("topright",names(large),lty=seq(1,ncol(large)),lwd=1,cex=0.75,
       col=seq(1,ncol(large)),bg="white",title="Fraction")
mtext("C",outer=TRUE,adj=0,line=-29.4)

dev.off()
