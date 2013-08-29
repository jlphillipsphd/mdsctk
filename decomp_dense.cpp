//
// 
//                This source code is part of
// 
//                        M D S C T K
// 
//       Molecular Dynamics Spectral Clustering ToolKit
// 
//                        VERSION 1.2.0
// Written by Joshua L. Phillips.
// Copyright (c) 2013, Joshua L. Phillips.
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

  // Declare the supported options.
  po::options_description cmdline_options;
  po::options_description program_options("Program options");
  program_options.add_options()
    ("help,h", "show this help message and exit")
    ("apsp-file,f", po::value<string>(&apsp_filename)->default_value("apsp.dat"), "Input:  Lower triangle of dense matrix file (string:filename)")
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
  cout << "evals-file =     " << evals_filename << endl;
  cout << "evecs-file =     " << evecs_filename << endl;
  cout << endl;

  /* Locals */
  int n, lda, info, lwork;
  int n_sym = 0;
  double wkopt;
  double* work;
  /* Local arrays */
  char jobz[8] = "Vectors";
  char uplo[6] = "Lower";

  // File input streams
  ifstream distances;

  // File output streams
  ofstream eigenvalues;
  ofstream eigenvectors;

  // Open files...
  distances.open(apsp_filename.c_str());
  eigenvalues.open(evals_filename.c_str());
  eigenvectors.open(evecs_filename.c_str());

  // Determine size of the matrix in the file
  distances.seekg(0,ios::end);
  n_sym = (distances.tellg() * sizeof(char) / sizeof(double));
  distances.seekg(0,ios::beg);
  for (n = 3; (n*(n-1)/2)+n < n_sym; n++);
  lda = n;

  // Allocate data structures...
  double *w = new double[n];
  double *a = new double[n*lda];
  double g = 0.0;

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

  // Query and allocate the optimal workspace
  lwork = -1;
  dsyev_( jobz, uplo, &n, a, &lda, w, &wkopt, &lwork, &info );
  lwork = (int) wkopt;
  work = new double[lwork];

  // Eigen decomposition...
  dsyev_( jobz, uplo, &n, a, &lda, w, work, &lwork, &info );

  if ( info > 0 ) {
    cerr << "ERROR: The algorithm failed to converge for " << info << " eigenvalues." << endl;
    return -1;
  }
  else {
    cout << "Number of converged eigenvalues/vectors found: "
	 << n << endl;

  }

  for (int x = n-1; x >= 0; x--) {

    // Scale the vector
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

  eigenvalues.close();
  eigenvectors.close();

  delete [] work;
  delete [] w;
  delete [] a;

  return 0;
}
