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

void affinity(CSC_matrix &A, double sigma_a = 1.0) {
  cout << "Sigma: " << sigma_a << endl;
  // Make affinity matrix...
  for (int x = 0; x < A.n; x++)
    for (int y = A.pcol[x]; y < A.pcol[x+1]; y++)
      A.M[y] = exp(-(A.M[y] * A.M[y]) / (2.0 * sigma_a * sigma_a)); 
}

void affinity(CSC_matrix &A, int k_a, double K, bool pSet = false) {

  double *sigma_a = new double[A.n];

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
    for (int y = 0; y < k_a && y < sorted_A[x].size(); y++)
      sigma_a[x] += sorted_A[x][y];
    sigma_a[x] /= (double) k_a;
  }
  if (pSet) {
    entropic_affinity_sigmas(A.n, k_a, K, sorted_A, sigma_a);
  }
  delete [] sorted_A;

  // Make affinity matrix...
  for (int x = 0; x < A.n; x++)
    for (int y = A.pcol[x]; y < A.pcol[x+1]; y++)
      A.M[y] = exp(-(A.M[y] * A.M[y]) / (2.0 * sigma_a[x] * sigma_a[A.irow[y]])); 

  for (int x = 1; x < A.n; x++)
    sigma_a[0] += sigma_a[x];
  cout << "Average sigma: " << (sigma_a[0] / (double) A.n) << endl;
  cout << endl;

  delete [] sigma_a;

}

void normalize(CSC_matrix &A) {
  // Turn distances into normalized affinities...
  double *d_a = new double[A.n];

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
  
  delete [] d_a;
}

int main(int argc, char* argv[])
{

  const char* program_name = "auto_heir_decomp_sparse";
  bool optsOK = true;
  gmx::initForCommandLine(&argc,&argv);
  copyright(program_name);
  cout << "   Reads the symmetric CSC format sparse matrix from" << endl;
  cout << "   input-file, and heirarchically decomposes the " << endl;
  cout << "   Laplacian matrix until relaxation time convergence" << endl;
  cout << "   criteria are met as the following reference:" << endl;
  cout << "   [1] B. Nadler and M. Galun, \"Fundamental Limitations" << endl;
  cout << "   of Spectral Clustering,\" in Advances in Neural Information" << endl;
  cout << "   Processing Systems 19, 2007, pp. 1017–1024." << endl;
  cout << "   eigenvalues/vectors of the normalized laplacian" << endl;
  cout << endl;
  cout << "   Use -h or --help to see the complete list of options." << endl;
  cout << endl;

  // Option vars...
  int k_a;
  double sigma;
  int nev = 2;
  double c1 = 1.2;
  double c2 = 2.0;
  double K;
  bool pSet = false;
  string ssm_filename;
  string output_filename;
  string ndx_filename;
  string residuals_filename;

  // Declare the supported options.
  po::options_description cmdline_options;
  po::options_description program_options("Program options");
  program_options.add_options()
    ("help,h", "show this help message and exit")
    ("sigma,s", po::value<double>(&sigma)->default_value(1.0), "Input:  Kernel sigma (double)")
    ("relaxation,r", po::value<double>(&c1)->default_value(1.2), "Input:  Relaxation cutoff parameter, c1 (double)")
    ("partition,p", po::value<double>(&c2)->default_value(2.0), "Input:  Partition cutoff parameter, c2 (double)")
    ("ssm-file,f", po::value<string>(&ssm_filename)->default_value("distances.ssm"), "Input:  Symmetric sparse matrix file (string:filename)")
    ("output,o", po::value<string>(&output_filename)->default_value("clusters.dat"), "Output:  Cluster assignment file (string:filename)")
    ("ndx,n", po::value<string>(&ndx_filename)->default_value("clusters.ndx"), "Output: Cluster assignment index file (string:filename)")    
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
  cout << "sigma =      " << sigma << endl;
  cout << "ssm-file =   " << ssm_filename << endl;
  cout << "output =     " << output_filename << endl;
  cout << "ndx =        " << ndx_filename << endl;
  cout << endl;

  // Stacks
  vector<vector<int> > work;
  vector<vector<int> > completed;

  // Defining variables;
  double  *Ax;  // Array for residual calculation
  double  residual = 0.0;
  double  max_residual = 0.0;

  // SSM Matrix
  CSC_matrix A(ssm_filename);
  Ax = new double[A.n];

  // File output streams
  ofstream output;
  ofstream ndx;

  // EPS
  double eps = getEPS();

  // Open files
  output.open(output_filename.c_str());
  ndx.open(ndx_filename.c_str());

  // Get affinities
  affinity(A,sigma);

  // Setup work
  work.resize(1);
  work[0].resize(A.n);
  for (int x = 0; x < A.n; x++)
    work[0][x] = x;

  while (work.size()) {
    vector<int> current = work[work.size()-1];
    work.pop_back();
    CSC_matrix current_A;
    A.syslice(current,current_A);
    normalize(current_A);

    double* d = NULL; // values
    double* Z = NULL; // vectors
    double nevm = runARPACK(nev,current_A,d,Z); 
    cout << "Number of converged eigenvalues/vectors found: "
	 <<  nevm << endl;
  
    int *labels = new int[current_A.n];
    kmeans(current_A.n,nev,nev,Z,labels);
    // kmeans(current_A.n,1,nev,&Z[A.n],labels);
    
    // Get slice indices
    vector<int> islice1 = select(current,0,labels);
    vector<int> islice2 = select(current,1,labels);

    if (islice1.size() == 0 || islice2.size() == 0) {
      cout << "Defunct partition..." << endl;
    }

    // Get slices
    CSC_matrix slice1;
    A.syslice(islice1,slice1);
    normalize(slice1);
    CSC_matrix slice2;
    A.syslice(islice2,slice2);
    normalize(slice2);
    
    double *slice1_d = NULL;
    double *slice1_Z = NULL;
    int nev1 = runARPACK(nev,slice1,slice1_d,slice1_Z);
    cout << "Number of converged eigenvalues/vectors found: "
	 <<  nev1 << endl;
    
    double *slice2_d = NULL;
    double *slice2_Z = NULL;
    int nev2 = runARPACK(nev,slice2,slice2_d,slice2_Z);
    cout << "Number of converged eigenvalues/vectors found: "
	 << nev2 << endl;

    double current_t,slice1_t,slice2_t,ratio;
    current_t = (1.0/(1.0-d[0]));
    slice1_t = (1.0/(1.0-slice1_d[0]));
    slice2_t = (1.0/(1.0-slice2_d[0]));
    ratio = slice1_t / slice2_t;
    if (ratio < 1.0)
      ratio = slice2_t / slice1_t;

    cout << "Main:   " << current_t << endl;
    cout << "Slice1: " << slice1_t << endl;
    cout << "Slice2: " << slice2_t << endl;
    cout << "Relaxation: " << c1*(slice1_t + slice2_t) << endl;
    cout << "Partition:  " << ratio << endl;

    if (nev1+nev2+nevm!=6) {
      cout << "No convergence. Skipping decomposition..." << endl;
      completed.push_back(current);
    }
    else {
      if (current_t < c1*(slice1_t + slice2_t)) {
	completed.push_back(current);
      }
      else if (ratio > c2) {
	if (slice1_t > slice2_t) {
	  work.push_back(islice1);
	  completed.push_back(islice2);
	}
	else {
	  work.push_back(islice2);
	  completed.push_back(islice1);
	}
      }
      else {
	work.push_back(islice1);
	work.push_back(islice2);
      }
    }

    delete [] labels;
    delete [] slice1_d;
    delete [] slice1_Z;
    delete [] slice2_d;
    delete [] slice2_Z;
    delete [] d;
    delete [] Z;
  }

  int *clusters = new int[A.n];
  cout << "Number of clusters: " << completed.size() << endl;
  cout << endl;

  for (int x = 0; x < completed.size(); x++) {
    int idx = 0;
    ndx << "[cluster_" << x+1 << "]" << endl;
    for (int y = 0; y < completed[x].size(); y++) {
      ndx << completed[x][y]+1 << " ";
      clusters[completed[x][y]] = x+1;
      if (++idx > 19) {
	ndx << endl;
	idx = 0;
      }
    }
    ndx << endl;
    ndx << endl;
  }

  for (int x = 0; x < A.n; x++)
    output << clusters[x] << endl;

  ndx.close();
  output.close();

  delete [] clusters;
  delete [] Ax;
  return 0;

} // main.
