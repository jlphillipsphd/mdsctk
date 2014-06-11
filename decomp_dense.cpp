//
// 
//                This source code is part of
// 
//                        M D S C T K
// 
//       Molecular Dynamics Spectral Clustering ToolKit
// 
//                        VERSION 1.2.2
// Written by Joshua L. Phillips.
// Copyright (c) 2012-2014, Joshua L. Phillips.
// check out http://github.com/douradopalmares/mdsctk/ for more information.
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
// 
// If you want to redistribute modifications, please consider that
// derived work must not be called official MDSCTK. Details are found
// in the README & LICENSE files - if they are missing, get the
// official version at github.com/douradopalmares/mdsctk/.
// 
// To help us fund MDSCTK development, we humbly ask that you cite
// the papers on the package - you can find them in the top README file.
// 
// For more info, check our website at http://github.com/douradopalmares/mdsctk/
// 
//

// Local
#include "config.h"
#include "mdsctk.h"

#ifdef HAVE_MAGMA
#include <magma.h>
#endif

int main(int argc, char* argv[]) {

  const char* program_name = "decomp_dense";
  bool optsOK = true;
  copyright(program_name);
  cout << "   Reads the lower-triangle of a dense matrix from" << endl;
  cout << "   input-file, and computes the number of requested" << endl;
  cout << "   eigenvalues/vectors of the double-centered, squared" << endl;
  cout << "   distances using LAPACK (Metric MDS)." << endl;
  cout << endl;
  cout << "   Use -h or --help to see the complete list of options." << endl;
  cout << endl;

  // Option vars...
  string apsp_filename;
  string evals_filename;
  string evecs_filename;
  string residuals_filename;

  // Declare the supported options.
  po::options_description cmdline_options;
  po::options_description program_options("Program options");
  program_options.add_options()
    ("help,h", "show this help message and exit")
    ("apsp-file,f", po::value<string>(&apsp_filename)->default_value("apsp.dat"), "Input:  Lower triangle of dense matrix file (string:filename)")
    ("residuals-file,r", po::value<string>(&residuals_filename)->default_value("residuals.dat"), "Output: Residuals file (string:filename)")    
    ("evals-file,v", po::value<string>(&evals_filename)->default_value("eigenvalues.dat"), "Output:  Eigenvalues file (string:filename)")
    ("evecs-file,e", po::value<string>(&evecs_filename)->default_value("eigenvectors.dat"), "Output: Eigenvectors file (string:filename)")    
    ;
  cmdline_options.add(program_options);

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
  po::notify(vm);    

  if (vm.count("help")) {
    cout << "usage: " << program_name << " [options]" << endl;
    cout << cmdline_options << endl;
    return 1;
  }

  if (!optsOK) {
    return -1;
  }

  cout << "Running with the following options:" << endl;
  cout << "apsp-file =      " << apsp_filename << endl;
  cout << "residuals-file = " << residuals_filename << endl;
  cout << "evals-file =     " << evals_filename << endl;
  cout << "evecs-file =     " << evecs_filename << endl;
  cout << endl;

  /* Locals */
  int n, lda, info, lwork, liwork, iwkopt;
  int n_sym = 0;
  double wkopt;
  double* work;
  int* iwork;
  /* Local arrays */
  char jobz[8] = "Vectors";
  char uplo[6] = "Lower";
  char notrans = 'N';

  // File input streams
  ifstream distances;

  // File output streams
  ofstream eigenvalues;
  ofstream eigenvectors;
  ofstream residuals;

  // Open files...
  distances.open(apsp_filename.c_str());
  eigenvalues.open(evals_filename.c_str());
  eigenvectors.open(evecs_filename.c_str());
  residuals.open(residuals_filename.c_str());

  // Determine size of the matrix in the file
  distances.seekg(0,ios::end);
  n_sym = (distances.tellg() * sizeof(char) / sizeof(double));
  distances.seekg(0,ios::beg);
  for (n = 3; (n*(n-1)/2)+n < n_sym; n++);
  lda = n;

  // Allocate data structures...
  double *w = new double[n];
  double *a = new double[n*lda];
  double *z = new double[n*lda];
  double *ax = new double[n];
  double g = 0.0;
  double residual = 0.0;
  double max_residual = 0.0;
  double eps = getEPS();

  // Read in the symmetric matrix...
  for (int x = 0; x < n; x++) {
    distances.read((char*) &a[x*n+x], (sizeof(double) / sizeof(char)) * (n-x));
    for (int y = x; y < n; y++) {
      // Squared distances...
      a[(x*n)+y] *= a[(x*n)+y];
      w[x] += a[(x*n)+y];
      w[y] += a[(x*n)+y];
    }
  }
  distances.close();

  // Row/Column means
  for (int x = 0; x < n; x++) {
    w[x] /= n;
    g += w[x];
  }
  // Global mean
  g /= n;

  // Double centering and scaling
  for (int x = 0; x < n; x ++)
    for (int y = x; y < n; y++) {
      a[(x*n)+y] += g - w[x] - w[y];
      a[(x*n)+y] *= -0.5;
    }
  memcpy(z,a,sizeof(double)*n*lda);

  // Query and allocate the optimal workspace
  lwork = -1;
  liwork = -1;
#ifdef HAVE_MAGMA
  magma_init();
  magma_print_devices();
  magma_dsyevd( jobz[0], uplo[0], n, a, lda, w, &wkopt, lwork, &iwkopt, liwork, &info );
#else
  dsyev_( jobz, uplo, &n, a, &lda, w, &wkopt, &lwork, &info );
  iwkopt = 1.0;
#endif
  lwork = (int) wkopt;
  work = new double[lwork];
  liwork = (int) iwkopt;
  iwork = new int[liwork];

  // Eigen decomposition...
#ifdef HAVE_MAGMA
  magma_dsyevd( jobz[0], uplo[0], n, a, lda, w, work, lwork, iwork, liwork, &info );
#else
  dsyev_( jobz, uplo, &n, a, &lda, w, work, &lwork, &info );
#endif

  if ( info > 0 ) {
    cerr << "ERROR: The algorithm failed to converge for " << info << " eigenvalues." << endl;
    return -1;
  }
  else {
    cout << "Number of converged eigenvalues/vectors found: "
	 << n << endl;

  }

  for (int x = n-1; x >= 0; x--) {

    int i = 1;
    double t = 1.0;
    double b = 0.0;
#ifdef HAVE_MAGMA
    magma_dsymv(uplo[0], n, t,
		z, lda, &a[n*x], i,
		b, ax, i);
    cout << ax[0] << endl;
    t = -w[x];
    magma_daxpy(n, t, &a[n*x], i, ax, i);
    residual = magma_dnrm2(n,ax,i)/fabs(w[x]);
#else
    dsymv_(uplo, &n, &t, z,
	   &lda, &a[n*x], &i, &b,
	   ax, &i);
    t = -w[x];
    daxpy_(&n, &t, &a[n*x], &i, ax, &i);
    residual = dnrm2_(&n, ax, &i)/fabs(w[x]);
#endif

    if (residual > max_residual)
      max_residual = residual;
#ifdef DECOMP_WRITE_DOUBLE
    residuals.write((char*) &residual, sizeof(double) / sizeof(char));
#else
    residuals << residual << endl;
#endif

    // Scale the vectors
    for (int y = 0; y < n; y++)
      a[n*x+y] *= sqrt(fabs(w[x]));

#ifdef DECOMP_WRITE_DOUBLE
    eigenvalues.write((char*) &w[x],(sizeof(double) / sizeof(char)));
    eigenvectors.write((char*) &a[n*x],(sizeof(double) * n) / sizeof(char));
#else
    eigenvalues << w[x] << endl;
    for (int y = 0; y < n; y++)
      eigenvectors << a[(n*x)+y] << " ";
    eigenvectors << endl;
#endif
  }

  cout << "Maximum residual: " << max_residual
       << " (eps: " << eps << ")" << endl;
  if (max_residual > eps) {
    cout << "*** Max residual too high (max_r > eps)!" << endl;
    cout << "*** Please, check results manually..." << endl;
  }

  eigenvalues.close();
  eigenvectors.close();
  residuals.close();

  delete [] work;
  delete [] iwork;
  delete [] w;
  delete [] a;
  delete [] ax;
  delete [] z;

  return 0;
}
