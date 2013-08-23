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

1.4 GSL
    You will need to install the GNU Scientific Library.
    (Debian/Ubuntu Packages: libgsl0ldbl and libgsl0-dev)
    http://www.gnu.org/software/gsl/

1.5 BLAS/ATLAS
    You will need to install BLAS and/or ATLAS for your system.
    (Debian/Ubuntu Packages: libblas3gf and libblas-dev)
    http://www.netlib.org/blas/ or http://math-atlas.sourceforge.net/

1.6 Berkeley (C++) DB Library
    You will need to install libdb, libdb_cxx, and related
    C++ development headers.
    (Debian/Ubuntu Packages: libdb++-dev and dependencies)
    http://www.oracle.com/technetwork/products/berkeleydb/downloads/index.html
 
1.7 R
    You will need to install R to perform the final stages of the
    clustering process automatically, or produce replicate-cluster
    assignment plots (using the 'fields' package which is available on
    the CRAN after installation).
    (Debian/Ubuntu Packages: r-base and r-base-dev)
    http://www.r-project.org/

2 Compilation

2.1 Unpacking/Compilation
    $ tar xzf mdsctk-1.0.0.tar.gz
    $ cd mdsctk-1.0.0
    $ cmake .
    $ make

2.2 Setup environment (bash example follows...)
    $ export MDSCTK_HOME=path to MDSCTK
    $ export PATH=${PATH}:${MDSCTK_HOME}

All binaries should now be built. You will need to set the
environment variable MDSCTK_HOME to the main source directory. You
might also consider putting the source directory in your PATH for ease
of use. Setting MDSCTK_HOME is NECESSARY for some scripts/programs,
but having the source directory in your PATH is optional.

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

1. auto_decomp_sparse

   Performs self-tuning specral decomposition of the graph laplacian
   based off of ideas developed in:
   Zelnik-manor, L., & Perona, P. (2005). Self-tuning spectral
   clustering. Advances in Neural Information Processing Systems 17
   (Vol. 2, pp. 1601–1608). MIT Press. Retrieved from
   http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.84.7940

   Usage: ./auto_decomp_sparse [# eigenvalues/vectors] [k_sigma]
      Reads the symmetric CSC format sparse matrix from the files
      sym_distances.dat, sym_row_indices.dat, & sym_col_indices.dat
      and computes the number of requested eigenvalues/vectors of the
      normalized laplacian using ARPACK and a gaussian kernel of width
      sigma, where sigma is calculated using the average k-smallest
      values in each column.

2. auto_decomp_sparse_nystrom

   Performs the same computation as auto_decomp_sparse and takes the
   same command-line arguments, but also uses the Nystrom method to
   project out-of-sample points obtained from the general CSC format
   sparse matrix in nonsym_distances.dat, nonsym_row_indices.dat, and
   nonsym_col_indices.dat.
   
3. bb_xtc_to_phipsi

   Usage: ./bb_xtc_to_phipsi [xtc file]...  Convert the provided xtc
      file to phipsi angles and write the results to standard output.

4. check_xtc

   Usage: ./check_xtc [xtc file]
      Report stats on the provided xtc file.

5. clustering_histogram.r

   Reads the cluster assignment data from cluster.dat and the
   replicate assignment data from assignment.dat, and outputs the
   corresponding joint replicate-cluster probability histogram in
   histogram.dat.

6. clustering_nmi.r

   Reads the joint replicate-cluster probability histogram from
   histogram.dat and prints the normalized mutual information for the
   histogram.

7. decomp_sparse

   Performs standard specral decomposition of the graph laplacian as
   developed in:
   Weiss, Y. (1999). Segmentation using eigenvectors: a unifying
   view. Proceedings of the Seventh IEEE International Conference on
   Computer Vision (pp. 975–982). IEEE. doi:10.1109/ICCV.1999.790354

   Usage: ./decomp_sparse [# eigenvalues/vectors] [sigma]
      Reads the symmetric CSC format sparse matrix from the files
      sym_distances.dat, sym_row_indices.dat, & sym_col_indices.dat
      and computes the number of requested eigenvalues/vectors of the
      normalized laplacian using a gaussian kernal of width sigma and
      ARPACK.

8. decomp_sparse_nystrom

   Performs the same computation as decomp_sparse and takes the same
   command-line arguments, but also uses the Nystrom method to project
   out-of-sample points obtained from the general CSC format sparse
   matrix in nonsym_distances.dat, nonsym_row_indices.dat, and
   nonsym_col_indices.dat.


9. density.r

   Usage: density.r [k] <sigma>
      Computes the kernel density estimate for the sparse distances
      in distances.dat. The number of k nearest neighbors is required
      but the bandwidth of the kernel, sigma, can be supplied or
      guesstimated based the data.

10. entropy.r

   Usage: entropy.r [k]
      Computes the local entropy of the given sparse
      matrix with indices from indices.dat and the densities
      in density.dat. The number of nearest neighbors,
      k, is required.

11. kmeans.r

   Usage: kmeans.r [k]
      Performs standard k-means clustering on the provided
      eigenvectors from (eigenvalues.dat and eigenvectors.dat).
      The number of clusters requested (k), can be 2>=k<=nev
      where nev is the number of eigenvectors. The results are
      written to clusters.dat, and a breakdown of assignments
      by cluster is written to clusters.ndx.

12. knn_data

    Usage: ./knn_data [# threads] [k] [vector size] [fitting data
       file] Computes the k nearest neighbors of all pairs of vectors
       in the given binary data files.

13. knn_rms

    Usage: ./knn_rms [# threads] [k] [topology file] [fitting xtc file]
       Computes the k nearest neighbors of all pairs of structures in
       the given xtc file. A topology PDB file should be provided for
       determining the mass of each atom.

    Output is a matrix of sorted distances and matrix of corresponding
    indices for each distance.

14. make_sysparse

    A symmetric CSC matrix is constructed from the data in
    distances.dat and indices.dat. The result is placed in
    sym_distances.dat, sym_row_indices.dat, and sym_col_indices.dat.

    Usage: ./make_sysparse [k] <output k>
       Converts the results from knn_rms into CSC format.

    Normally, the number of nearest neighbors in the input distances
    is used for constructing the CSC matrix.  However, you can set
    <output k> <= [k] in order to subselect the number of neighbors to
    consider in the CSC representation. This makes it easy to store a
    large number of neighbors using knn_* but then use a subset for,
    say, computing approximate geodesic distances.

15. make_gesparse

    A general CSC matrix is constructed from the data in distances.dat
    and indices.dat. The result is placed in nonsym_distances.dat,
    nonsym_row_indices.dat, and nonsym_col_indices.dat.

    Usage: ./make_gesparse [k] <output k>
       Converts the results from knn_rms into CSC format.

    Normally, the number of nearest neighbors in the input distances
    is used for constructing the CSC matrix.  However, you can set
    <output k> <= [k] in order to subselect the number of neighbors to
    consider in the CSC representation. This makes it easy to store a
    large number of neighbors using knn_* but then use a subset for,
    say, computing approximate geodesic distances.

16. phipsi_to_sincos

    Converts the angles in phipsi.dat into polar coordinate
    representation (sincos.dat), which are appropriate vectors for
    distance calculations in torsion angle space.

17. plot_histogram.r

    Reads the cluster assignment data from cluster.dat and the
    replicate assignment data from assignment.dat, and outputs the
    corresponding joint replicate-cluster probability histogram as EPS
    file using the 'fields' package in R.

18. probability.r

    Usage: density.r [k] <sigma>
       Computes the kernel density estimate for the sparse distances
       in distances.dat, and convert the results to estimated
       probailities.
       The number of k nearest neighbors is required but the bandwidth of
       the kernel, sigma, can be supplied or guesstimated based the data.

19. rms_test

    Computes RMSD between the provided reference structure and all of
    the structures in the XTC file. This is a nice way to verify/test
    distance comparisons, but note that different codes produce
    slightly different results (usually within ~5% difference).

    Usage: ./rms_test [reference structure] [fitting xtc file]
       Computes the RMSD of all structures for the given in the xtc
       file. A template structure should be provided as the reference
       structure.

20. split_xtc

    A utility to sample from an XTC file, producing suitable files for
    using the Nystrom out-of-sample projection method.

    Usage: ./split_xtc [xtc file] [n]
       Splits the provided xtc file into two separate.  xtc files
       (landmarks.xtc and remainder.xtc) where landmarks.xtc will
       contain every n-th frame from the input xtc file, while
       remainder.xtc will contain all other frames from the input
       trajectory.

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

1. knn_rms [# threads] [k] [topology file] [fitting xtc file]

   Input: 2 100 trp-cage.pdb trp-cage.xtc
   Output: distances.dat indices.dat

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

2. make_sysparse [k] <output k>

   Input: 100 (distances.dat indices.dat)
   Output: sym_row_indices.dat sym_col_indices.dat sym_distances.dat

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

3. auto_decomp_sparse [# eigenvalues/vectors] [k]

   Input: 10 10 (sym_row_indices.dat sym_col_indices.dat
   sym_distances.dat)
   Output: eigenvalues.dat eigenvectors.dat residuals.dat

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

4. kmeans.r

   Input: (eigenvectors.dat)
   Output: clusters.dat

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

6. clustering_histogram.r

   Input: (clusters.dat assignment.dat)
   Output: histogram.dat

   This script takes the cluster assignment and frame assignment data
   and calculates the joint replicate-cluster assignment
   probablility, p(R,C). In other words, what fraction of the frames were
   assigned to cluster C and were also assigned to group R. This
   information tells us something about how trajectory is exploring
   the available conformation space. This information is analyzed in
   the following steps to make a quantitative prediction as to how
   well the different assignment groups overlap in conformation space.

   This histogram.dat file contains a CxR probability matrix.
 
7. clustering_nmi.r

   Input: (histogram.dat)
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

8. plot_histogram.r

   Input: (histogram.dat)
   Output: histogram.eps

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
