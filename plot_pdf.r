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
    program.name <- "plot_pdf.r"
    source(paste(Sys.getenv("MDSCTK_HOME"),"/mdsctk.r",sep=""))
}

cat("   Reads in a joint assignment-cluster probability distribution,\n")
cat("   plots the distribution as a heatmap, and also prints the\n")
cat("   normalized mutual information between assignments and clusters\n")
cat("   in the lower-right corner of the plot.\n")
cat("\n")
cat("   Use -h or --help to see the complete list of options.\n")
cat("\n")

parser$add_argument("-p","--pdf",default="pdf.dat",
                    help="Joint PDF file [default %(default)s]")
parser$add_argument("-o","--output",default="pdf.eps",
                    help="Joint PDF file [default %(default)s]")

myargs <- parser$parse_args()

cat("Running with the following options:\n")
cat(paste("pdf =    ",myargs$pdf,"\n"))
cat(paste("output = ",myargs$output,"\n"))
cat("\n")

if (!suppressPackageStartupMessages(require("fields",character.only=TRUE))) {
    cat("   Please install the R package 'fields' to use this MDSCTK R scripts.\n")
    cat("\n")
    q()
}

data <- as.matrix(read.table(myargs$pdf,header=FALSE))
n <- nrow(data)
m <- ncol(data)
mycol <- rev(heat.colors(50))
nmi <- normalmutualinf(data)
mark <- "A"

postscript(myargs$output,width=5.5,height=4.6,onefile=FALSE,horizontal=FALSE)

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
temp <- dev.off()
