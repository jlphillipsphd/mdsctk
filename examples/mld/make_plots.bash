#!/bin/bash
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
## official version at <website>.
## 
## To help us fund MDSCTK development, we humbly ask that you cite
## the papers on the package - you can find them in the top README file.
## 
## For more info, check our website at http://github.com/douradopalmares/mdsctk/
## 
##

WINDOW=500

./figure-k1.r semirigid-unsmoothed.dat semirigid-smoothed.dat figure-03.eps

./figure-k2.r halffolded-unsmoothed.dat halffolded-smoothed.dat figure-04.eps

./figure-k1.r correlated-unsmoothed.dat correlated-smoothed.dat figure-05.eps

./figure-pw1.r correlated-unsmoothed-small.dat correlated-unsmoothed-medium.dat \
    correlated-unsmoothed-large.dat ${WINDOW} figure-06.eps

./figure-pw2.r correlated-smoothed-small.dat correlated-smoothed-medium.dat \
    correlated-smoothed-large.dat ${WINDOW} figure-07.eps

./figure-prot.r gb1-rmsd.dat gb1-unsmoothed-estimates.dat \
    gb1-smoothed-estimates.dat figure-11.eps

./figure-prot.r trp-cage-rmsd.dat trp-cage-unsmoothed-estimates.dat \
    trp-cage-smoothed-estimates.dat figure-12.eps

./figure-prot.r nsp1-rmsd.dat nsp1-unsmoothed-estimates.dat \
    nsp1-smoothed-estimates.dat figure-13.eps

./figure-prot.r nup116-rmsd.dat nup116-unsmoothed-estimates.dat \
    nup116-smoothed-estimates.dat figure-14.eps
