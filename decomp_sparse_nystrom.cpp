//
// 
//                This source code is part of
// 
//                        M D S C T K
// 
//       Molecular Dynamics Spectral Clustering ToolKit
// 
//                        VERSION 1.1.1
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

// Local
#include "config.h"

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

void sp_dgemv(int n, int *irow, int *pcol, double *A,
	      double *v, double *w) {
  
  int i,j,k;
  
  for (i=0; i<n; i++) w[i] = 0.0;
  
  for (i=0; i<n; i++) {
    k = pcol[i];
    for (j=k; j<pcol[i+1]; j++) {
      w[i] += v[irow[j]]*A[j];
    }
  }
}

int main(int argc, char* argv[])
{

  if (argc != 3) {
    cerr << endl;
    cerr << "   MDSCTK " << MDSCTK_VERSION_MAJOR << "." << MDSCTK_VERSION_MINOR << endl;
    cerr << endl;
    cerr << "Usage: " << argv[0] << " [# eigenvalues/vectors] [sigma]" << endl;
    cerr << "   Reads the symmetric CSC format sparse matrix from the" << endl;
    cerr << "   files sym_distances.dat, sym_row_indices.dat, &" << endl;
    cerr << "   sym_col_indices.dat and computes the number of" << endl;
    cerr << "   requested eigenvalues/vectors of the normalized" << endl;
    cerr << "   laplacian using a gaussian kernal of width sigma" << endl;
    cerr << "   and ARPACK. Also, the nonsymmetric CSC format sparse" << endl;
    cerr << "   matrix in nonsym_distances.dat, nonsym_row_indices.dat, &" << endl;
    cerr << "   nonsym_col_indices.dat is projected onto the earlier" << endl;
    cerr << "   eigen vectors (i.e. out-of-sample prediction)." << endl;
    cerr << endl;
    return -1;
  }

  int     nev = atoi(argv[1]); // Number of eigenvectors...
  double  sigma_a = atof(argv[2]); // Sigma
  
  // General
  int     n;   // Dimension of the problem.
  int     m;   // Outer dimension

  // Main affinity matrix
  int     nnzA;
  int     *irowA;
  int     *pcolA;
  double  *A;   // Pointer to an array that stores the lower
		// triangular elements of A.

  // Expanded affinity matrix
  int     nnzB;
  int     *irowB;
  int     *pcolB;
  double  *B;   // Pointer to an array that stores the
		// sparse elements of B.

  // File input streams
  ifstream distances;
  ifstream row_pointers;
  ifstream col_pointers;

  // File output streams
  ofstream eigenvalues;
  ofstream eigenvectors;
  ofstream residuals;

  // EPS
  double eps = 1.0;
  do { eps /= 2.0; } while (1.0 + (eps / 2.0) != 1.0);
  eps = sqrt(eps);
 
  // Open files
  distances.open("sym_distances.dat");
  row_pointers.open("sym_row_indices.dat");
  col_pointers.open("sym_col_indices.dat");
  eigenvalues.open("eigenvalues.dat");
  eigenvectors.open("eigenvectors.dat");
  residuals.open("residuals.dat");

  // Determine size of col pointer file
  col_pointers.seekg(0,ios::end);
  n = (col_pointers.tellg() * sizeof(char) / sizeof(int));
  col_pointers.seekg(0,ios::beg);
  pcolA = new int[n];
  col_pointers.read((char*) pcolA, (sizeof(int) / sizeof(char)) * n);
  n--;
  nnzA = pcolA[n];
  A = new double[nnzA];
  irowA = new int[nnzA];
  row_pointers.read((char*) irowA, (sizeof(int) / sizeof(char)) * nnzA);
  distances.read((char*) A, (sizeof(double) / sizeof(char)) * nnzA);

  distances.close();
  row_pointers.close();
  col_pointers.close();

  // Open files
  distances.open("nonsym_distances.dat");
  row_pointers.open("nonsym_row_indices.dat");
  col_pointers.open("nonsym_col_indices.dat");

  col_pointers.seekg(0,ios::end);
  m = (col_pointers.tellg() * sizeof(char) / sizeof(int));
  col_pointers.seekg(0,ios::beg);
  pcolB = new int[m];
  col_pointers.read((char*) pcolB, (sizeof(int) / sizeof(char)) * m);
  m--;
  nnzB = pcolB[m];
  B = new double[nnzB];
  irowB = new int[nnzB];
  row_pointers.read((char*) irowB, (sizeof(int) / sizeof(char)) * nnzB);
  distances.read((char*) B, (sizeof(double) / sizeof(char)) * nnzB);

  distances.close();
  row_pointers.close();
  col_pointers.close();

  // Turn distances into normalized affinities...
  double *d_a = new double[n];
  double *d_b = new double[m];

  // Make affinity matrices...
  for (int x = 0; x < nnzA; x++)
    A[x] = exp(-(A[x] * A[x]) / (2.0 * sigma_a * sigma_a));
  
  for (int x = 0; x < nnzB; x++)
    B[x] = exp(-(B[x] * B[x]) / (2.0 * sigma_a * sigma_a));
  
  // Calculate D_A
  for (int x = 0; x < n; x++)
    d_a[x] = 0.0;
  for (int x = 0; x < n; x++) {
    for (int y = pcolA[x]; y < pcolA[x+1]; y++) {
      d_a[x] += A[y];
      d_a[irowA[y]] += A[y];
    }
  }
  for (int x = 0; x < n; x++)
    d_a[x] = 1.0 / sqrt(d_a[x]);

  // Calculate D_B
  for (int x = 0; x < m; x++)
    d_b[x] = 0.0;
  for (int x = 0; x < m; x++) {
    for (int y = pcolB[x]; y < pcolB[x+1]; y++) {
      d_b[x] += B[y];
    }
  }
  for (int x = 0; x < m; x++)
    d_b[x] = 1.0 / sqrt(d_b[x]);

  // Normalize the affinity matrix...
  for (int x = 0; x < n; x++) {
    for (int y = pcolA[x]; y < pcolA[x+1]; y++) {
      A[y] *= d_a[irowA[y]] * d_a[x];
    }
  }

  // Normalized B matrix...
  for (int x = 0; x < m; x++) {
    for (int y = pcolB[x]; y < pcolB[x+1]; y++) {
      B[y] *= d_a[irowB[y]] * d_b[x];
    }
  }

  delete [] d_a;
  delete [] d_b;


  // Eigen decomposition of nomalized affinity matrix...

  // ARPACK setup...
  double  *Ax = new double[n];  // Array for residual calculation
  double  residual = 0.0;
  double  max_residual = 0.0;
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
  double *extrap_evec = new double[m];
  double *norm_evec = new double[n];

  while (ido != 99) {
    dsaupd_(&ido, &bmat, &n, which,
	    &nev, &tol, resid,
	    &ncv, &V[1], &ldv,
	    &iparam[1], &ipntr[1], &workd[1],
	    &workl[1], &lworkl, &info);
    
    if (ido == -1 || ido == 1) {
      // Matrix-vector multiplication
      sp_dsymv(n, irowA, pcolA, A,
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
#endif

    // Extrapolate remaining points onto the vector space
    for (int y = 0; y < n; y++)
      norm_evec[y] = Z[(n*x)+y] / d[x];
    sp_dgemv(m, irowB, pcolB, B, norm_evec, extrap_evec);

#ifdef DECOMP_WRITE_DOUBLE
    eigenvectors.write((char*) extrap_evec,(sizeof(double) * m) / sizeof(char));
#else
    for (int y = 0; y < m; y++)
      eigenvectors << extrap_evec[y] << " ";
    eigenvectors << endl;
#endif

    // Calculate residual...
    // Matrix-vector multiplication
    sp_dsymv(n, irowA, pcolA, A,
	     &Z[n*x], Ax);

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

  cout << "Max residual: " << max_residual
       << " (eps: " << eps << ")" << endl;
  if (max_residual > eps) {
    cout << "*** Sum of residuals too high (max_r > eps)!" << endl;
    cout << "*** Please, check results manually..." << endl;
  }

  eigenvalues.close();
  eigenvectors.close();
  residuals.close();

  delete [] irowA;
  delete [] pcolA;
  delete [] A;

  delete [] irowB;
  delete [] pcolB;
  delete [] B;

  // ARPACK
  delete [] lselect;
  delete [] d;
  delete [] resid;
  delete [] Ax;
  delete [] V;
  delete [] iparam;
  delete [] ipntr;
  delete [] workd;
  delete [] workl;
  delete [] norm_evec;
  delete [] extrap_evec;

  return 0;
}
