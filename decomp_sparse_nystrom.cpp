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

int main(int argc, char* argv[])
{

  const char* program_name = "decomp_sparse_nystrom";
  bool optsOK = true;
  copyright(program_name);
  cout << "   Reads the symmetric CSC format sparse matrix from" << endl;
  cout << "   input-file, and computes the number of requested" << endl;
  cout << "   eigenvalues/vectors of the normalized laplacian" << endl;
  cout << "   using ARPACK and a gaussian kernel of width sigma." << endl;
  cout << "   The general CSC format sparse matrix is projected" << endl;
  cout << "   onto the eigenvectors of the symmetric matrix for" << endl;
  cerr << "   out-of-sample prediction." << endl;
  cout << endl;
  cout << "   Use -h or --help to see the complete list of options." << endl;
  cout << endl;

  // Option vars...
  double sigma_a;
  int nev;
  string ssm_filename;
  string gsm_filename;
  string evals_filename;
  string evecs_filename;
  string residuals_filename;

  // Declare the supported options.
  po::options_description cmdline_options;
  po::options_description program_options("Program options");
  program_options.add_options()
    ("help,h", "show this help message and exit")
    ("sigma,q", po::value<double>(&sigma_a), "Input:  Standard deviation of gaussian kernel (real)")
    ("nevals,n", po::value<int>(&nev), "Input:  Number of eigenvalues/vectors (int)")
    ("ssm-file,s", po::value<string>(&ssm_filename)->default_value("distances.ssm"), "Input:  Symmetric sparse matrix file (string:filename)")
    ("gsm-file,g", po::value<string>(&gsm_filename)->default_value("distances.gsm"), "Input:  General sparse matrix file (string:filename)")
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
  if (!vm.count("sigma")) {
    cout << "ERROR: --sigma not supplied." << endl;
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
  cout << "sigma =          " << sigma_a << endl;
  cout << "nevals =         " << nev << endl;
  cout << "ssm-file =       " << ssm_filename << endl;
  cout << "gsm-file =       " << gsm_filename << endl;
  cout << "evals-file =     " << evals_filename << endl;
  cout << "evecs-file =     " << evecs_filename << endl;
  cout << "residuals-file = " << residuals_filename << endl;
  cout << endl;

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
  ifstream ssm;
  ifstream gsm;

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
  gsm.open(gsm_filename.c_str());
  eigenvalues.open(evals_filename.c_str());
  eigenvectors.open(evecs_filename.c_str());
  residuals.open(residuals_filename.c_str());

  // Read symmetric CSC matrix
  ssm.read((char*) &n, (sizeof(int) / sizeof(char)));
  pcolA = new int[n+1];
  ssm.read((char*) pcolA, (sizeof(int) / sizeof(char)) * (n+1));
  nnzA = pcolA[n];
  A = new double[nnzA];
  irowA = new int[nnzA];
  ssm.read((char*) irowA, (sizeof(int) / sizeof(char)) * nnzA);
  ssm.read((char*) A, (sizeof(double) / sizeof(char)) * nnzA);
  ssm.close();

  // Read general CSC matrix
  gsm.read((char*) &m, (sizeof(int) / sizeof(char)));
  pcolB = new int[m+1];
  gsm.read((char*) pcolB, (sizeof(int) / sizeof(char)) * (m+1));
  nnzB = pcolB[m];
  B = new double[nnzB];
  irowB = new int[nnzB];
  gsm.read((char*) irowB, (sizeof(int) / sizeof(char)) * nnzB);
  gsm.read((char*) B, (sizeof(double) / sizeof(char)) * nnzB);
  gsm.close();

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
