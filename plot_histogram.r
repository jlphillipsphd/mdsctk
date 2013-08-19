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

require(fields)

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

data <- as.matrix(read.table("histogram.dat",header=FALSE))
n <- nrow(data)
m <- ncol(data)
mycol <- rev(heat.colors(50))
nmi <- normalmutualinf(data)
mark <- "A"

postscript("histogram.eps",width=5.5,height=4.6,onefile=FALSE,horizontal=FALSE)

par(mar=c(3,3,0.5,4.6)+0.1,mgp=c(2.1,0.8,0))
image(seq(1,n),seq(1,m),data,zlim=c(0,max(data)),col=mycol,
      ylab="Cluster",xlab="Replicate")
box()
# mtext(mark,outer=TRUE,side=3,adj=0,line=-1.5,cex=2)
mtext(sprintf("NMI=%4.3f",nmi),outer=TRUE,side=1,adj=1,line=-1,cex=2)
image.plot(seq(1,n),seq(1,m),data,zlim=c(0,1.0/max(c(n,m))),col=mycol,
           legend.only=TRUE,
           legend.args=list(cex=0.75,text=expression(p(R,C))),
           legend.mar=4.1,useRaster=TRUE)
dev.off()
