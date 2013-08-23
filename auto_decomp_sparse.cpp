//
// 
//                This source code is part of
// 
//                        M D S C T K
// 
//       Molecular Dynamics Spectral Clustering ToolKit
// 
//                        VERSION 1.1.2
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

// Standard
// C
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
// C++
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

// Boost
#include <boost/program_options.hpp>

// Local
#include "config.h"
#include "mdsctk.h"

// ARPACK routines
extern "C" {
  
  void dsaupd_(int *ido, char *bmat, int *n, char *which,
	       int *nev, double *tol, double *resid,
	       int *ncv, double *V, int *ldv,
	       int *iparam, int *ipntr, double *workd,
	       double *workl, int *lworkl, int *info);
  
  void dseupd_(int *rvec, char *HowMny, int *select,
	       double *d, double *Z, int *ldz,
	       double *sigma, char *bmat, int *n,
	       char *which, int *nev, double *tol,
	       double *resid, int *ncv, double *V,
	       int *ldv, int *iparam, int *ipntr,
	       double *workd, double *workl,
	       int *lworkl, int *info);

  // For checking residuals
  void daxpy_(const int *n, const double *da, const double *dx,
	      const int *incx, double *dy, const int *incy);
  double dnrm2_(const int *n, const double *dx, const int *incx);
  
} // end FORTRAN definitions

namespace po = boost::program_options;
using namespace std;

// Sparse Routines
void sp_dsymv(int n, int *irow, int *pcol, double *A,
	      double *v, double *w) {

  int i,j,k;
  double t = 0.0;
  
  for (i=0; i<n; i++) w[i] = 0.0;

  for (i=0; i<n; i++) {
    t = v[i];
    k = pcol[i];
    if ((k!=pcol[i+1])&&(irow[k]==i)) {
      w[i] += t*A[k];
      k++;
    }
    for (j=k; j<pcol[i+1]; j++) {
      w[irow[j]] += t*A[j];
      w[i] += v[irow[j]]*A[j];
    }
  }
}

int main(int argc, char* argv[])
{

  const char* program_name = "auto_decomp_parse";
  bool optsOK = true;
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
    cout << endl;
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

  if (!optsOK) {
    return -1;
  }

  cout << "Running with the following options:" << endl;
  cout << "k-sigma =    " << k_a << endl;
  cout << "nevals =     " << nev << endl;
  cout << "ssm-file =   " << ssm_filename << endl;
  cout << "evals-file = " << evals_filename << endl;
  cout << "evecs-file = " << evecs_filename << endl;
  cout << endl;

  // Defining variables;
  int     n;   // Dimension of the problem.
  int     nnz;
  int     *irow;
  int     *pcol;
  double  *A;   // Pointer to an array that stores the lower
		// triangular elements of A.
  double  *Ax;  // Array for residual calculation
  double  residual = 0.0;
  double  max_residual = 0.0;

  // File input streams
  ifstream ssm;

  // File output streams
  ofstream eigenvalues;
  ofstream eigenvectors;
  ofstream residuals;

  // EPS
  double eps = 1.0;
  do { eps /= 2.0; } while (1.0 + (eps / 2.0) != 1.0);
  eps = sqrt(eps);

  // Open files
  ssm.open(ssm_filename.c_str());
  eigenvalues.open(evals_filename.c_str());
  eigenvectors.open(evecs_filename.c_str());
  residuals.open(residuals_filename.c_str());

  // Read symmetric CSC matrix
  ssm.read((char*) &n, (sizeof(int) / sizeof(char)));
  pcol = new int[n+1];
  ssm.read((char*) pcol, (sizeof(int) / sizeof(char)) * (n+1));
  nnz = pcol[n];
  Ax = new double[n];
  A = new double[nnz];
  irow = new int[nnz];
  ssm.read((char*) irow, (sizeof(int) / sizeof(char)) * nnz);
  ssm.read((char*) A, (sizeof(double) / sizeof(char)) * nnz);
  ssm.close();

  // Begin AFFINTY

  // Turn distances into normalized affinities...
  double *sigma_a = new double[n];
  double *d_a = new double[n];

  // Calculate sigmas...
  vector<double> *sorted_A = new vector<double>[n];
  for (int x = 0; x < n; x++)
    sorted_A[x].clear();
  for (int x = 0; x < n; x++)
    for (int y = pcol[x]; y < pcol[x+1]; y++) {
      sorted_A[x].push_back(A[y]);
      sorted_A[irow[y]].push_back(A[y]);
    }
  for (int x = 0; x < n; x++) {
    partial_sort(sorted_A[x].begin(),
		 sorted_A[x].begin() +
		 ((sorted_A[x].size() < k_a)?sorted_A[x].size():k_a),
		 sorted_A[x].end());
    sigma_a[x] = 0;
    for (int y = 0; y < k_a && y < sorted_A[x].size(); y++)
      sigma_a[x] += sorted_A[x][y];
    sigma_a[x] /= (double) k_a;
  }
  delete [] sorted_A;

  // Make affinity matrix...
  for (int x = 0; x < n; x++)
    for (int y = pcol[x]; y < pcol[x+1]; y++)
      A[y] = exp(-(A[y] * A[y]) / (2.0 * sigma_a[x] * sigma_a[irow[y]]));

  // Calculate D_A
  for (int x = 0; x < n; x++)
    d_a[x] = 0.0;
  for (int x = 0; x < n; x++) {
    for (int y = pcol[x]; y < pcol[x+1]; y++) {
      d_a[x] += A[y];
      d_a[irow[y]] += A[y];
    }
  }
  for (int x = 0; x < n; x++)
    d_a[x] = 1.0 / sqrt(d_a[x]);

  // Normalize the affinity matrix...
  for (int x = 0; x < n; x++) {
    for (int y = pcol[x]; y < pcol[x+1]; y++) {
      A[y] *= d_a[irow[y]] * d_a[x];
    }
  }

  delete [] sigma_a;
  delete [] d_a;

  // End AFFINITY

  // ARPACK variables...
  int ido = 0;
  char bmat = 'I';
  char which[2];
  which[0] = 'L';
  which[1] = 'A';
  double tol = 0.0;
  double *resid = new double[n];
  // NOTE: Need about one order of magnitude more arnoldi vectors to
  // converge for the normalized Laplacian (according to residuals...)
  int ncv = ((10*nev+1)>n)?n:(10*nev+1);
  double *V = new double[(ncv*n)+1];
  int ldv = n;
  int *iparam = new int[12];
  iparam[1] = 1;
  iparam[3] = 100 * nev;
  iparam[4] = 1;
  iparam[7] = 1;
  int *ipntr = new int[15];
  double *workd = new double[(3*n)+1];
  int lworkl = ncv*(ncv+9);
  double *workl = new double[lworkl+1];
  int info = 0;
  int rvec = 1;
  char HowMny = 'A';
  int *lselect = new int[ncv];
  double *d = new double[nev];
  double *Z = &V[1];
  int ldz = n;
  double sigma = 0.0;
 
  while (ido != 99) {
    dsaupd_(&ido, &bmat, &n, which,
	    &nev, &tol, resid,
	    &ncv, &V[1], &ldv,
	    &iparam[1], &ipntr[1], &workd[1],
	    &workl[1], &lworkl, &info);
    
    if (ido == -1 || ido == 1) {
      // Matrix-vector multiplication
      sp_dsymv(n,irow,pcol,A,
	       &workd[ipntr[1]],
	       &workd[ipntr[2]]);
    }
  }
    
  dseupd_(&rvec, &HowMny, lselect,
	  d, Z, &ldz,
	  &sigma, &bmat, &n,
	  which, &nev, &tol,
	  resid, &ncv, &V[1],
	  &ldv, &iparam[1], &ipntr[1],
	  &workd[1], &workl[1],
	  &lworkl, &info);

  cout << "Number of converged eigenvalues/vectors found: "
       << iparam[5] << endl;
  
  for (int x = nev-1; x >= 0; x--) {
#ifdef DECOMP_WRITE_DOUBLE
    eigenvalues.write((char*) &d[x],(sizeof(double) / sizeof(char)));
    eigenvectors.write((char*) &Z[n*x],(sizeof(double) * n) / sizeof(char));
#else
    eigenvalues << d[x] << endl;
    for (int y = 0; y < n; y++)
      eigenvectors << Z[(n*x)+y] << " ";
    eigenvectors << endl;
#endif

    // Calculate residual...
    // Matrix-vector multiplication
    sp_dsymv(n,irow,pcol,A,
	     &Z[n*x],
	     Ax);

    double t = -d[x];
    int i = 1;
    daxpy_(&n, &t, &Z[n*x], &i, Ax, &i);
    residual = dnrm2_(&n, Ax, &i)/fabs(d[x]);
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

  delete [] irow;
  delete [] pcol;
  delete [] Ax;
  delete [] A;

  // ARPACK
  delete [] lselect;
  delete [] d;
  delete [] resid;
  delete [] V;
  delete [] iparam;
  delete [] ipntr;
  delete [] workd;
  delete [] workl;

  return 0;

} // main.
