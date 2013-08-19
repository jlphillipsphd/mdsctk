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

valin <- "eigenvalues.dat"
vecin <- "eigenvectors.dat"
resin <- "residuals.dat"

e.values <- as.double(scan(valin,quiet=TRUE))
e.vectors <- matrix(scan(vecin,quiet=TRUE),ncol=length(e.values))
e.residuals <- as.double(scan(resin,quiet=TRUE))

set.seed(0) # Change for different results...
clusters <- kmeans(e.vectors,
                   ncol(e.vectors),
                   iter.max=30,
                   nstart=10)$cluster

myout <- file("clusters.dat","w")
write(clusters,myout,ncolumns=1)
close(myout)

myout <- file("clusters.ndx","w")
for (n in seq(1,ncol(e.vectors))) {
  writeLines(sprintf("[cluster_%d]",n),con=myout)
  write(which(clusters==n),myout,ncolumns=20)
  writeLines("",con=myout)
}
close(myout)

q(save="no")
