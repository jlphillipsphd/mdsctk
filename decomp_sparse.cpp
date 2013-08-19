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
    cerr << "   and ARPACK." << endl;
    cerr << endl;
    return -1;
  }

  // Defining variables;
  int     nev = atoi(argv[1]); // Number of eigenvectors...
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
  pcol = new int[n];
  col_pointers.read((char*) pcol, (sizeof(int) / sizeof(char)) * n);
  n--;
  nnz = pcol[n];
  Ax = new double[n];
  A = new double[nnz];
  irow = new int[nnz];
  row_pointers.read((char*) irow, (sizeof(int) / sizeof(char)) * nnz);
  distances.read((char*) A, (sizeof(double) / sizeof(char)) * nnz);

  distances.close();
  row_pointers.close();
  col_pointers.close();

  // Begin AFFINTY

  // Turn distances into normalized affinities...
  double sigma_a = atof(argv[2]);
  double *d_a = new double[n];

  // ofstream naff;
  // naff.open("sym_n_affinities.dat");

  // Make affinity matrix...
  for (int x = 0; x < nnz; x++)
    A[x] = exp(-(A[x] * A[x]) / (2.0 * sigma_a * sigma_a));

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

  cout << "Max residual: " << max_residual
       << " (eps: " << eps << ")" << endl;
  if (max_residual > eps) {
    cout << "*** Sum of residuals too high (max_r > eps)!" << endl;
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
