mdsctk
======

MDSCTK - Molecular Dynamics Spectral Clustering Toolkit

Version: 1.2.0
Date: 2013-08-20

Author:
Joshua L. Phillips <jphillips@lanl.gov>
T-6/CNLS, Los Alamos National Laboratory

Thanks for considering MDSCTK for your spectral analysis needs!

The development of MDSCTK is mainly funded by academic research grants.
To help us fund development, we humbly ask that you cite the MDSCTK papers:

* Validating clustering of molecular dynamics simulations using
  polymer models
  J. L. Phillips, M. E. Colvin and S. Newsam
  BMC Bioinformatics, 12(1), 445. (2011)
  doi:10.1186/1471-2105-12-445

* Analyzing dynamical simulations of intrinsically disordered proteins
  using spectral clustering
  J. L. Phillips, M. E. Colvin, E. Y. Lau and S. Newsam
  Proc. of the 2008 IEEE Int. Conf. on Bioinf. and Biomed. Workshops
  (pp. 17–24). Philadelphia, PA: IEEE. (2008)
  doi:10.1109/BIBMW.2008.4686204

*************
GENERAL INFO:
*************

MDSCTK is free software, distributed under the GNU General Public License. 

If you want to distribute a modified version or use part of MDSCTK
in your own program, remember that the entire modified code must be licensed 
under GPL, and that it must clearly be labeled as derived work. It should 
not use the name "MDSCTK", and make sure support questions are
directed to you instead of the MDSCTK developers.

The MDSCTK is a suite of tools for performing clustering or
dimensionality reduction of molecular dynamics simulations using
spectral methods. The toolkit consists of a set of programs and
scripts that work in a pipelined fashion, where the output from one
program or script becomes input for another. This approach allows the
user to perform processing on intermediate data or just use the parts
of the pipeline which are useful for a particular task, and ignore the
rest.

Here are some of the things that MDSCTK does well:

1. Parallel computation of RMSD or vector distance matrices.

2. Sparse eigen decomposition on Compressed Sparse Column matrices.

3. Simple Phi-Psi angle extraction from backbone coordinates.

4. Nystrom approximation technique to cluster out-of-sample data.

New features:

1. Basic implementation of Isomap.

2. Sparse contact map calculation.

3. Knn codes for sparse vectors (complements the contact maps).

4. New codes for estimating entropy, probability, and density.

5. Overall more consistent interface among codes/scripts.

6. Bash environment setup script.

************
INSTALLATION
************

1 Prerequisites

1.1 CMake
    You will need cmake to generate a Makefile in order to complete
    compilation/installation.
    (Debian/Ubuntu Packages: cmake and dependecies)
    http://www.cmake.org/

1.2 GROMACS
    You will need to make sure that you have a working installation of
    GROMACS and that you have sourced the GMXRC file to set up your
    GROMACS environment variables properly.
    (Debian/Ubuntu Packages: gromacs and gromacs-dev)
    http://www.gromacs.org/

1.3 ARPACK
    You will need to install the ARPACK sparse eigen solver library.
    (Debian/Ubuntu Packages: libarpack2 and libarpack2-dev)
    http://www.caam.rice.edu/software/ARPACK/

1.4 BLAS/ATLAS
    You will need to install BLAS and/or ATLAS for your system.
    (Debian/Ubuntu Packages: libblas3gf and libblas-dev)
    http://www.netlib.org/blas/ or http://math-atlas.sourceforge.net/

1.5 Berkeley (C++) DB Library
    You will need to install libdb, libdb_cxx, and related
    C++ development headers.
    (Debian/Ubuntu Packages: libdb++-dev and dependencies)
    http://www.oracle.com/technetwork/products/berkeleydb/downloads/index.html

1.6 Boost C++ Graph and Program Options Libraries
    You will need to install libboost-program-options-dev and
    libboost-graph-dev at least, but I would recommend just installing
    libboost-dev which will grab all of the boost development files.
    (Debian/Ubuntu Packages: libboost-program-options-dev,
    libboost-graph-dev, and dependencies (or just install the entire
    boost development files for simplicity: libboost-dev)
    http://www.boost.org/
 
1.7 R
    You will need to install R to perform the final stages of the
    clustering process automatically, or produce replicate-cluster
    assignment plots (using the 'fields' and 'argparse' packages which
    are available on the CRAN after installation).
    (Debian/Ubuntu Packages: r-base and r-base-dev)
    http://www.r-project.org/

2 Compilation

2.1 Unpacking/Compilation
    $ tar xzf mdsctk-1.0.0.tar.gz
    $ cd mdsctk-1.0.0
    $ cmake .
    $ make

2.2 Setup environment (bash example follows...)
    $ source MDSCTK.bash

You will need to set the environment variable MDSCTK_HOME to the main
source directory manually if you are not using bash. You might also
consider putting the source directory in your PATH for ease of
use. Setting MDSCTK_HOME is NECESSARY for some scripts/programs, but
having the source directory in your PATH is optional.

There are a lot of uni/linux configurations out there, and I am not
familiar with all of them. If you experience problems, more than
likely it is because a prerequisite isn't installed properly or it is
installed in an non-standard location. You will need to modify some
CMAKE variables using ccmake or cmake-gui to help the compiler find
the needed library/header files.

(If you gave difficulties with setting up the build system using
cmake, an auxiliary Makefile (Makefile.old) is provided with the
source code and can be edited by hand. Sometimes, this will be easier
than trying to get cmake to find all of the required elements using
it's search tools. However, the provided Makefile is currently
deprecated, and will eventually be removed from the project. Just copy
Makefile.old to Makefile and also config.h.old to config.h, then make
any needed edits to Makefile.)

*******************
Basic Documentation
*******************

All binaries and R scripts accept the -h and --help options at the
command line, which will provide basic usage information. Because of
this, provided here is simply a list of the available tools with any
additional comments that are not found in the program descriptions.

1. auto_decomp_sparse

   Performs self-tuning specral decomposition of the graph laplacian
   based off of ideas developed in:
   Zelnik-manor, L., & Perona, P. (2005). Self-tuning spectral
   clustering. Advances in Neural Information Processing Systems 17
   (Vol. 2, pp. 1601–1608). MIT Press. Retrieved from
   http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.84.7940
   or alternatively for entropic affinities:
   M. Vladymyrov and M. A. Carreira-Perpiñán, “Entropic Affinities:
   Properties and Efficient Numerical Computation,” in Proceedings of
   the 30th International Conference on Machine Learning (ICML-13),
   2013, vol. 28, no. 3, pp. 477–485. Retrieved from
   http://jmlr.org/proceedings/papers/v28/vladymyrov13.pdf

2. auto_decomp_sparse_nystrom
   
3. bb_xtc_to_phipsi

4. check_xtc

5. clustering_nmi.r

6. clustering_pdf.r

7. contact_profile

8. decomp_dense

   Performs standard metric scaling of a dense matrix.

9. decomp_sparse

   Performs standard specral decomposition of the graph laplacian as
   developed in:
   Weiss, Y. (1999). Segmentation using eigenvectors: a unifying
   view. Proceedings of the Seventh IEEE International Conference on
   Computer Vision (pp. 975–982). IEEE. doi:10.1109/ICCV.1999.790354

10. decomp_sparse_nystrom

11. density.r

12. dijkstra

    Computes all pairs of shortest paths for a sparse CSC matrix,
    yeilding a dense matrix for performing metric scaling using
    decomp_dense. This is the ISOMAP algorithm as developed in:
    J. B. Tenenbaum, V. de Silva, and J. C. Langford, “A global
    geometric framework for nonlinear dimensionality reduction,”
    Science, vol. 290, no. 5500, pp. 2319–2323, Dec. 2000.

13. entropy.r

14. kmeans.r

15. knn_data

16. knn_data_ocl

    Proof-of-concept CPU/GPU acceleration using OpenCL.

17. knn_data_sparse

18. knn_rms

19. make_sysparse

20. make_gesparse

21. phipsi_to_sincos

22. plot_pdf.r

23. probability.r

24. rms_test

25. split_xtc

********
EXAMPLES
********

The examples/ directory contains scripts that show how to perform
various tasks using MDSCTK. Some familiarity with GROMACS analysis
tools is assumed. If you are not familiar with GROMACS, then you might
want to read the GROMACS documentation on how one uses the trjconv
tool to convert your trajectories into GROMACS XTC format.

To test your installation, just try:
$ cd examples/
$ ./cluster_rms.bash

This script performs autoscaled spectral clustering on the test data
found in examples/trp-cage.pdb and examples/trp-cage.xtc

Each command in the script performs a specific task that creates the
input files for the next command. Follow along by opening
cluster_rms.bash in your favorite editor:

1. knn_rms -t [# threads] -k [k] -p [topology file] -r [fitting xtc file] -d [distances file] -i [indices file]

   Input: -t 2 -k 100 -p trp-cage.pdb -r trp-cage.xtc
   Output: -d distances.dat -i indices.dat (defaults)

   This program performs parallel computation of the RMS distances
   between all pairs of structures in the provided XTC file.

   The number of threads used in the computation can be specified,
   which might not scale very well past 8 cores in general since disk
   I/O becomes a larger bottleneck than RMS computations. It is only
   set to 2 in the script, but feel free to experiment. I usually get
   max throughput using 8 cores on most modern systems.

   We are only keeping the k-nearest neighbors to keep data storage at
   a minimum. Picking k is difficult, but it should be much larger
   than the number of neighbors used for computing scaling factors
   (see auto_decomp_sparse below). In general, k will be more-or-less
   constrained by the amount of disk space needed to store the
   matrix. If we have N frames in the trajectory, then we will have to
   store N*k double values on-disk, and the corresponding integer
   indices (N*k integer values). However, numbers around 50-100 should
   work well for most problems (just a rule of thumb, so please try a
   few values!).

   The topology information is extracted from the provided PDB file (any
   GROMACS-compatible structure file should work though) to obtain the
   mass of each atom for the weighted RMSD fits. The coordinates are
   not used.

   The XTC file should contain all of the structures that you want to
   cluster.

   When this command completes, you will have two new files in the
   current directory: distances.dat and indices.dat. These files
   contain the k-nearest RMS distances for each frame, and the
   corresponding indices of those k-nearest frames,
   respectively. These are raw double and integer files, respectively,
   if one is interested in using the output with other
   programs/scripts.

2. make_sysparse -k [k] -n [output k] -d [distance file] -i [index file] -o [CSC matrix file]

   Input: -k 100 (distances.dat indices.dat are defaults)
   Output: distances.ssm

   This program converts the data files from the last step
   (distances.dat and indices.dat) into symmetric Compressed Sparse
   Column format, which is the input format for the next program in the
   pipeline.

   Again, the number k should be the same as the last step so that
   make_sparse will produce the correct CSC matrix geometry. However,
   you can provide a value for <output k> which is less than the one
   used above if want to try a subset of the stored values. This is
   useful for testing if your clustering is too sensitve to using a
   few less neighbors. In that case, it may be wise to increase k in
   the previous knn_* calculations until dropping a few neighbors
   causes no stability/sensitivity issues.	

   For details on the CSC format, refer to:
   http://netlib.org/linalg/html_templates/node92.html#SECTION00931200000000000000

3. auto_decomp_sparse -n [# eigenvalues/vectors] -k [k-sigma] -s [CSC matrix file] -v [eigenvalues file] -e [eigenvectors file] -r [residuals file]

   Input: -n 10 -k 10 (distance.ssm by default)
   Output: eigenvalues.dat eigenvectors.dat residuals.dat (defaults)

   This program reads the CSC distance matrix from the sym_*.dat
   files and computes the graph laplacian using an automatically tuned
   gaussian kernel based on the supplied k=10. This should be much less
   than than the number of neighbors in the matrix (chosen above) to
   ensure that the values of the laplacian will slowly decay to
   zero.

   The largest 10 eigenvectors of the laplacian are
   computed and stored in eigenvectors.dat. Corresponding eigenvalues
   and residuals are stored in eigenvalues.dat and residuals.dat,
   respectively. One should check that the residuals are small
   (eps<1e-8). Large residuals could indicate that your chosen
   k-parameter values above were too extreme (small or large). These
   files are ASCII text so they should be easy to import into other
   programs as well.

   Critically, the NUMBER of eigenvalues/vectors corresponds to the
   number of CLUSTERS you want to partition your trajectory into. If
   you would like to test your system across a wide range of clusters,
   then just select the maximum that you would like to investigate and
   then modify the kmeans.r script (used below) to use just the first
   two eigenvectors for two clusters, three eigenvectors for three
   clusters, etc. The key idea is that you only need to perform this
   step once to get a sufficient number of vectors because they will
   NEVER change (unless you change the scaling parameter). However,
   keeping too many vectors will slow down the computation, and result
   in a large eigenvectors.dat file, so be reasonable.

4. kmeans.r -k [#clusters] -v [eigenvalues file] -e [eigenvectors file] -o [clusters file] -n [clusters index file]

   Input: -k 10 (eigenvalues.dat eigenvectors.dat by default)
   Output: clusters.dat clusters.ndx (again by default)

   This is an R script that performs k-means clustering on the
   resulting eigenvectors, and outputs the cluster assignments for
   each frame (one per line) in clusters.dat

   If one doesn't prefer to use R (Rscript), then a custom k-means
   program or other general k-means package could be used in place of
   this script.

   If one is only interested in the cluster assignment for each frame,
   then look no futher! The next steps perform some interesting
   analysis from the clustering results, but clusters.dat contains the
   cluster assignment for each 

5. Rscript

   Output: assignment.dat

   This command creates a replicate assignment vector and places it in
   the file assignment.dat. The idea here is that we need to tell the
   next set of scripts how the frames are "naturally" partitioned.

   In this example, the trajectory is a contiguous 1000 frames from
   a single MD trajectory (folding trp-cage), but we decide to
   partition it into 10 consecutive time intervals. Hence, the first
   100 frames are assigned to interval 1, the next 200 are assigned to
   interval 2, etc.

   Other ways are possible, such as 5 intervals of 400 frames each, so
   one can experiment. If your trajectory file consists of ten
   independent simulations which have been concatenated, then it might
   make sense to make the assignment based on this information. You
   can even just assign every frame to the same interval (1) by just
   creating an assignment.dat file that has N lines, each containing a
   1. This wouldn't be very interesting, but I hope that it is clear
   that this assignment is based off of some other practical reason
   for why certain frames may belong together.

6. clustering_pdf.r -a [assignment file] -c [clusters file] -o [pdf file]

   Input: (clusters.dat assignment.dat by default)
   Output: pdf.dat (again, by default)

   This script takes the cluster assignment and frame assignment data
   and calculates the joint replicate-cluster assignment
   probablility, p(R,C). In other words, what fraction of the frames were
   assigned to cluster C and were also assigned to group R. This
   information tells us something about how trajectory is exploring
   the available conformation space. This information is analyzed in
   the following steps to make a quantitative prediction as to how
   well the different assignment groups overlap in conformation space.

   This pdf.dat file contains a CxR probability matrix.
 
7. clustering_nmi.r -p [pdf file]

   Input: (pdf.dat by default)
   Output: Normalized Mutual Information - NMI

   This script takes the probablities above and computes the
   normalized mutual information of the replicate-cluster
   assignment. A value of 0 indicates that the system is completely
   mixed (frames from any one replicate are equally distrbuted across
   all clusters), while a value of 1 indicates that the system is NOT
   mixed at all (frames from any one replicate were placed in a unique
   cluster). Values in-between 0 and 1 quantify how strongly mixed or
   unmixed the conformations are. Therefore, NMI of p(R,C) can be used
   to assess simlation convergence, divergence, or just conformational
   evolution over time (as in the example).

8. plot_pdf.r -p [pdf file] -o [pdf eps]

   Input: (pdf.dat by default)
   Output: pdf.eps (again. by default)

   If you have the 'fields' package installed in R, then this script
   will create a 3D plot of p(R,C). This allows visual examination of
   the result.

   The result for the trp-cage simulation is clear, the system gets
   closer and closer to the folded state over time as indicated by
   darker values (higher probability) along the main diagonal. It
   finally folds in the last four intervals (and corresponding
   clusters 7,8,9,10), where the p(R,C) plot is more mixed. However,
   it also quickly  unfolds back to conformations seen earlier
   (cluster 4, regions 4,5) in the last interval (10). By increasing
   the number of assignment intervals or number of clusters, the
   granularity of the results can be improved, but general trends my
   begin to disappear. Exploration is needed to determine the best set
   of parameters to use in this case.

Additional considerations...

Another thing to keep in mind is that all molecules should be
processed for periodic boundary conditions, and made whole
before attempting to use them with MDSCTK.

There is also an example script for using polar-transformed phi-psi
angle Euclidean distance for the distance metric instead of RMSD
(cluster_phipsi.bash). The only thing to keep in mind here is that
only the backbone N-CA-C atoms should be retained in the trajectory if
this is desired. The trjconv tool in GROMACS makes it easy to extract
this subset of atoms.

There is a script, cluster_data.bash, which illustrates how to use the
toolkit to cluster arbitrary data vectors. In this example, it uses
two concentric circles of data points, and partitions them into two
clusters.

The cluster_*_nystrom.bash scripts illustrate how to perform the
Nystrom approximation on out-of-sample data points. This is useful for
VERY large data sets where it is intractable to solve the clustering
problem for all points directly. Instead, a subsample of the data set
large enough to capture all of the necessary structural classes is
used, and then the remaining data can be projected onto the
eigenvectors of this subsample laplacian, allowing cluster
classification of the entire data set. The workflow is similar to the
one used in the example above, except that the neighbor calculations
for the out-of-sample data are needed, and the appropriate MDSCTK
utilities (*_nystrom) are used.

The cluster_contacts*.bash scripts show how to produce and use sparse
input vectors (in this case, contact maps).

The isomap_data.bash script illustrates how to perform Isomap on a
data set, and can be produced from RMSD, Phi-Psi angle, or contact map
data as well. However, the Nystrom extension, and even basic landmark
MDS, has not been added to this feature just yet, so you will only be
able to run this algorithm on relatively small data sets.
