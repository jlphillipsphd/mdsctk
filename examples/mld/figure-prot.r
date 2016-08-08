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
## check out http://github.com/jlphillipsphd/mdsctk/ for more information.
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
## For more info, check our website at http://github.com/jlphillipsphd/mdsctk/
## 
##

myargs <- commandArgs(TRUE)
rmsd.file <- myargs[1]
unsmoothed.file <- myargs[2]
smoothed.file <- myargs[3]
figure.file <- myargs[4]

temperature <- read.table("temperature.txt",header=FALSE)
rmsd <- read.table(rmsd.file,header=FALSE)
unsmoothed <- scan(unsmoothed.file)
smoothed <- scan(smoothed.file)

## Clean-up data
unsmoothed <- cbind(rmsd[,1],unsmoothed)
smoothed <- cbind(rmsd[,1],smoothed)

## Make plot
postscript(figure.file,width=3.35,height=5.8,onefile=FALSE,horizontal=FALSE)
layout(matrix(c(1,2,3,4)),heights=c(1.5,1.5,1.5,1.8))

## First plot
par(mar=c(0,3.2,0.2,0)+0.1,mgp=c(1.5,0.5,0))
plot(temperature[seq(1,nrow(temperature),100),1]/1000.0,
     temperature[seq(1,nrow(temperature),100),2],
     type='l',xaxt='n',ylab="Temperature (K)",xlab="Time (ns)",
     cex.axis=0.75,cex.lab=1.2)
mtext("A",outer=TRUE,adj=0,side=3,line=-1.2)

## Second plot
par(mar=c(0,3.2,0,0)+0.1,mgp=c(1.5,0.5,0))
plot(rmsd[,1]/1000.0,runmed(rmsd[,2]/10,k=101),
     col="red",type='l',lty=1,lwd=1,
     xaxt='n',cex.axis=0.75,ylab=expression(paste("RMSD (",ring(A),")")),
     ylim=c(0,1),cex.lab=1.2)
mtext("B",outer=TRUE,adj=0,side=3,line=-11.8)

## Third plot
par(mar=c(0,3.2,0,0)+0.1,mgp=c(1.5,0.5,0))
plot(unsmoothed[,1]/1000.0,runmed(unsmoothed[,2],k=501),
     col="red",type='l',lty=1,lwd=1,
     xaxt='n',cex.axis=0.75,ylab="Estimated Dim.",
     ylim=c(0,40),cex.lab=1.2)
mtext("C",outer=TRUE,adj=0,side=3,line=-22.2)

## Fourth plot
par(mar=c(2.4,3.2,0,0)+0.1,mgp=c(1.5,0.5,0))
plot(smoothed[,1]/1000.0,runmed(smoothed[,2],k=501),
     col="red",type='l',lty=1,lwd=1,
     cex.axis=0.75,ylab="Estimated Dim.",ylim=c(0,40),
     xlab="Time (ns)",cex.lab=1.2)
mtext("D",outer=TRUE,adj=0,side=3,line=-32.6)

dev.off()
