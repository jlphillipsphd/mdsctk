//
// 
//                This source code is part of
// 
//                        M D S C T K
// 
//       Molecular Dynamics Spectral Clustering ToolKit
// 
//                        VERSION 1.2.5
//
// Written by Joshua L. Phillips.
// Copyright (c) 2012-2015, Joshua L. Phillips.
// Check out http://www.cs.mtsu.edu/~jphillips/software.html for more
// information.
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of the
// License, or (at your option) any later version.
// 
// If you want to redistribute modifications, please consider that
// derived work must not be called official MDSCTK. Details are found
// in the README & LICENSE files - if they are missing, get the
// official version at github.com/douradopalmares/mdsctk/.
// 
// To help us fund MDSCTK development, we humbly ask that you cite the
// papers on the package - you can find them in the top README file.
// 
// For more info, check our website at
// http://www.cs.mtsu.edu/~jphillips/software.html
// 
//

// Local
#include "config.h"
#include "mdsctk.h"

int main(int argc, char* argv[])
{

  const char* program_name = "auto_decomp_sparse";
  bool optsOK = true;
  gmx::initForCommandLine(&argc,&argv);
  copyright(program_name);
  cout << "   Reads the symmetric CSC format sparse matrix from" << endl;
  cout << "   input-file, and computes the number of requested" << endl;
  cout << "   eigenvalues/vectors of the normalized laplacian" << endl;
  cout << "   using ARPACK and a gaussian kernel of width sigma," << endl;
  cout << "   where sigma is calculated using the average of the" << endl;
  cout << "   k-closest neighbor distances in each column." << endl;
  cout << "   Note that k is normally smaller than the number of" << endl;
  cout << "   neighbors used to construct the sparse CSC matrix." << endl;
  cout << endl;
  cout << "   Use -h or --help to see the complete list of options." << endl;
  cout << endl;

  // Option vars...
  int k_a;
  int nev;
  double K;
  bool pSet = false;
  string ssm_filename;
  string evals_filename;
  string evecs_filename;
  string residuals_filename;

  // Declare the supported options.
  po::options_description cmdline_options;
  po::options_description program_options("Program options");
  program_options.add_options()
    ("help,h", "show this help message and exit")
    ("k-sigma,k", po::value<int>(&k_a), "Input:  K-nn to average for sigmas (int)")
    ("k-perplexity,K", po::value<double>(&K), "Input:  Desired perplexity within knn (real)")
    ("nevals,n", po::value<int>(&nev), "Input:  Number of eigenvalues/vectors (int)")
    ("ssm-file,s", po::value<string>(&ssm_filename)->default_value("distances.ssm"), "Input:  Symmetric sparse matrix file (string:filename)")
    ("evals-file,v", po::value<string>(&evals_filename)->default_value("eigenvalues.dat"), "Output:  Eigenvalues file (string:filename)")
    ("evecs-file,e", po::value<string>(&evecs_filename)->default_value("eigenvectors.dat"), "Output: Eigenvectors file (string:filename)")    
    ("residuals-file,r", po::value<string>(&residuals_filename)->default_value("residuals.dat"), "Output: Residuals file (string:filename)")    
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
  if (!vm.count("k-sigma")) {
    cout << "ERROR: --k-sigma not supplied." << endl;
    cout << endl;
    optsOK = false;
  }
  if (!vm.count("nevals")) {
    cout << "ERROR: --nevals not supplied." << endl;
    cout << endl;
    optsOK = false;
  }
  if (vm.count("k-perplexity")) {
    pSet = true;
    if (ceil(K) > (double) k_a) {
      cout << "WARNING: --k-perplexity (" << K
	   << ") is higher than --k-sigma (" << k_a
	   << ")" << endl;
    }
  }

  if (!optsOK) {
    return -1;
  }

  cout << "Running with the following options:" << endl;
  cout << "k-sigma =        " << k_a << endl;
  if (pSet)
    cout << "k-perplexity =   " << K << endl;
  cout << "nevals =         " << nev << endl;
  cout << "ssm-file =       " << ssm_filename << endl;
  cout << "residuals-file = " << residuals_filename << endl;
  cout << "evals-file =     " << evals_filename << endl;
  cout << "evecs-file =     " << evecs_filename << endl;
  cout << endl;

  // Defining variables;
  double  *Ax;  // Array for residual calculation
  double  residual = 0.0;
  double  max_residual = 0.0;

  // SSM Matrix
  CSC_matrix A(ssm_filename);
  Ax = new double[A.n];

  // File output streams
  ofstream eigenvalues;
  ofstream eigenvectors;
  ofstream residuals;

  // EPS
  double eps = getEPS();

  // Open files
  eigenvalues.open(evals_filename.c_str());
  eigenvectors.open(evecs_filename.c_str());
  residuals.open(residuals_filename.c_str());

  // Begin AFFINTY

  // Turn distances into normalized affinities...
  double *sigma_a = new double[A.n];
  double *d_a = new double[A.n];

  // Calculate sigmas...
  vector<double> *sorted_A = new vector<double>[A.n];
  for (int x = 0; x < A.n; x++)
    sorted_A[x].clear();
  for (int x = 0; x < A.n; x++)
    for (int y = A.pcol[x]; y < A.pcol[x+1]; y++) {
      sorted_A[x].push_back(A.M[y]);
      sorted_A[A.irow[y]].push_back(A.M[y]);
    }
  for (int x = 0; x < A.n; x++) {
    while (sorted_A[x].size() > k_a)
      sorted_A[x].pop_back();
    sort(sorted_A[x].begin(),sorted_A[x].end());
    sigma_a[x] = 0;
    for (int y = 0; y < sorted_A[x].size(); y++)
      sigma_a[x] += sorted_A[x][y];
    sigma_a[x] /= (double) k_a;
  }
  if (pSet)
    entropic_affinity_sigmas(A.n, k_a, K, sorted_A, sigma_a);
  delete [] sorted_A;

  // Make affinity matrix...
  for (int x = 0; x < A.n; x++)
    for (int y = A.pcol[x]; y < A.pcol[x+1]; y++)
      A.M[y] = exp(-(A.M[y] * A.M[y]) / (2.0 * sigma_a[x] * sigma_a[A.irow[y]]));

  // Calculate D_A
  for (int x = 0; x < A.n; x++)
    d_a[x] = 0.0;
  for (int x = 0; x < A.n; x++) {
    for (int y = A.pcol[x]; y < A.pcol[x+1]; y++) {
      d_a[x] += A.M[y];
      d_a[A.irow[y]] += A.M[y];
    }
  }
  for (int x = 0; x < A.n; x++)
    d_a[x] = 1.0 / sqrt(d_a[x]);

  // Normalize the affinity matrix...
  for (int x = 0; x < A.n; x++) {
    for (int y = A.pcol[x]; y < A.pcol[x+1]; y++) {
      A.M[y] *= d_a[A.irow[y]] * d_a[x];
    }
  }

  for (int x = 1; x < A.n; x++)
    sigma_a[0] += sigma_a[x];
  cout << "Average sigma: " << (sigma_a[0] / (double) A.n) << endl;
  cout << endl;

  delete [] sigma_a;
  delete [] d_a;

  // End AFFINITY

  double* d; // values
  double* Z; // vectors
  cout << "Number of converged eigenvalues/vectors found: "
       << runARPACK(nev,A,d,Z) << endl;
  
  for (int x = nev-1; x >= 0; x--) {
#ifdef DECOMP_WRITE_DOUBLE
    eigenvalues.write((char*) &d[x],(sizeof(double) / sizeof(char)));
    eigenvectors.write((char*) &Z[A.n*x],(sizeof(double) * A.n) / sizeof(char));
#else
    eigenvalues << d[x] << endl;
    for (int y = 0; y < A.n; y++)
      eigenvectors << Z[(A.n*x)+y] << " ";
    eigenvectors << endl;
#endif

    // Calculate residual...
    // Matrix-vector multiplication
    sp_dsymv(A.n,A.irow,A.pcol,A.M,
	     &Z[A.n*x],
	     Ax);

    double t = -d[x];
    int i = 1;
    daxpy_(&A.n, &t, &Z[A.n*x], &i, Ax, &i);
    residual = dnrm2_(&A.n, Ax, &i)/fabs(d[x]);
    if (residual > max_residual)
      max_residual = residual;
#ifdef DECOMP_WRITE_DOUBLE
    residuals.write((char*) &residual, sizeof(double) / sizeof(char));
#else
    residuals << residual << endl;
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


  // ARPACK
  delete [] Ax;
  delete [] d;
  delete [] Z;

  return 0;

} // main.
