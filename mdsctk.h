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
// check out http://cnls.lanl.gov/~jphillips/ for more information.
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
// 
// If you want to redistribute modifications, please consider that
// derived work must not be called official MDSCTK. Details are found
// in the README & COPYING files - if they are missing, get the
// official version at cnls.lanl.gov/~jphillips/.
// 
// To help us fund MDSCTK development, we humbly ask that you cite
// the papers on the package - you can find them in the top README file.
// 
// For more info, check our website at http://cnls.lanl.gov/~jphillips/
// 
//

#ifndef MDSCTK_H
#define MDSCTK_H

// Standard
// C
#include <stddef.h>
//#include <stdlib.h>
#include <math.h>
// C++
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <cstdlib>

// OpenMP
#include <omp.h>

// GROMACS
#include <gromacs/tpxio.h>
#include <gromacs/xtcio.h>
#include <gromacs/index.h>
#include <gromacs/do_fit.h>
#include <gromacs/centerofmass.h>

// Boost
#include <boost/program_options.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

// Berkeley DB
#include <db_cxx.h>

#define STRIFY(str) #str
#define STRINGIFY(str) STRIFY(str)
#define asize(a) (sizeof(a)/sizeof((a)[0]))
#define SQR(A) (A*A)
#define AVG(A,B) ((A+B)/2.0)
#define INDEX(SIZE,I,J) ((SIZE * J) + I)

#define BLKSIZE 65536

typedef ::real (*coord_array)[3];
namespace po = boost::program_options;
using namespace std;

extern const double RAD2DEG;
extern const size_t update_interval;

// Edges struct, comparison, and sorting routines
struct edge {
  int from;
  int to;  
};

// CSC Matrix
struct CSC_matrix {
  CSC_matrix();
  CSC_matrix(const string filename);
  CSC_matrix(const CSC_matrix&);
  ~CSC_matrix();
  CSC_matrix& operator=(const CSC_matrix&);
  int     n;   // Dimension of the matrix (#cols)
  int     nnz;
  int     *irow;
  int     *pcol;
  double  *M;   // Pointer to an array that stores the
		// elements of the matrix.
  //void geslice(vector<int>& r, CSC_matrix& csc);
  void syslice(vector<int>& rc, CSC_matrix& csc);
  double& operator[](int);
  void init();
  void cleanup();
  void copy(const CSC_matrix&);
};

// TOP_file
class TOP_file {
public:
  TOP_file(const string init_filename);
  ~TOP_file();

  int get_natoms();
  ::real* get_mass();
  coord_array get_frame_ptr();
  void center(coord_array frame);
  ::real rmsd(coord_array ref_frame, coord_array fit_frame);
  void com(coord_array frame, int n, atom_id index[], rvec com);
  void get_index(const string ndx_filename, int &ndx_n, int* &ndx_index, char* &ndx_group);
  
private:
  string filename;
  int natoms;
  ::real *mass;
  coord_array frame;
  t_topology top;
  int ePBC;
  matrix box;
  char buf[256];

  void read_topology();
};

// XTC file
class XTC_file {
public:
  XTC_file(const string init_filename);
  ~XTC_file();
  
  int get_natoms();
  float get_time();
  float get_prec();
  int get_step();
  coord_array get_next_frame();
  coord_array get_next_frame_ptr();

private:
  string filename;
  t_fileio *file;
  int step;
  float time;
  matrix box;
  float prec;
  gmx_bool bOK;
  int natoms;
  coord_array frame;
};

// Permutation template
template <class T> struct permutation {
  vector<T> data;
  vector<int> indices;
  void sort(int k = 0) {
    if (!k) {
      indices.resize(data.size());
      for (int x = 0; x < data.size(); x++)
	indices[x] = x;
      std::sort(indices.begin(),indices.end(),*this);
      std::sort(data.begin(),data.end());
    }
    else {
      indices.resize(data.size());
      for (int x = 0; x < data.size(); x++)
	indices[x] = x;
      std::partial_sort(indices.begin(),indices.begin()+k,
			indices.end(),*this);
      std::partial_sort(data.begin(),data.begin()+k,
			data.end());
    }
  }
  bool operator()(int left, int right) const { return data[left]<data[right]; }
};

void copyright(const char* program_name = NULL);

double getEPS();

// Sparse Routines
void sp_dsymv(int n, int *irow, int *pcol, double *A,
	      double *v, double *w);
void sp_dgemv(int n, int *irow, int *pcol, double *A,
	      double *v, double *w);

// Distance metrics
double euclidean_distance(int size, double* reference, double* fitting);

double correlation_distance(int size, double* reference, double* fitting);

double euclidean_distance_sparse(int ref_size, int* ref_index, double* ref_data,
				 int fit_size, int* fit_index, double* fit_data);

// Sigma calculation routine
void entropic_affinity_sigmas(int n, int k, double K,
			      std::vector<double>* A, double* s);

// Split edges
void split_edges(int current_index, Dbc *cursor,
		 vector<int> &indices, vector<double> &distances);
int compare_edge(Db *db, const Dbt *key1, const Dbt *key2);

// Math
::real theta(::real pos1[], ::real pos2[], ::real pos3[],
	     bool degrees);

void crossprod(::real C[],
	       ::real x1, ::real y1, ::real z1,
	       ::real x2, ::real y2, ::real z2);

::real torsion(::real pos1[], ::real pos2[], ::real pos3[], ::real pos4[], bool degrees);

void sample(int n, int k, int *sample);

void kmeans(int n, int d, int k, double *data, int *labels,
	    int nstarts = 25, int maxit = 30);

vector<int> select(vector<int> &set, int k, int *labels);

// Runtime
int runARPACK(int nev, CSC_matrix &A, double* &d, double* &Z);
int runARPACK2(int nev, CSC_matrix &A, double* &d, double* &Z);

// ARPACK, BLAS, LAPACK, etc.
extern "C" {

  // ARPACK
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

  void dnaupd_(int *ido, char *bmat, int *n, char *which,
	       int *nev, double *tol, double *resid,
	       int *ncv, double *V, int *ldv,
	       int *iparam, int *ipntr, double *workd,
	       double *workl, int *lworkl, int *info);

  void dneupd_(int *rvec, char *HowMny, int *select,
	       double *dr, double* di, double *Z, int *ldz,
	       double *sigmar, double *sigmai, double* workev,
	       char *bmat, int *n,
	       char *which, int *nev, double *tol,
	       double *resid, int *ncv, double *V,
	       int *ldv, int *iparam, int *ipntr,
	       double *workd, double *workl,
	       int *lworkl, int *info);

  // LAPACK
  void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda,
	      double* w, double* work, int* lwork, int* info );

  // BLAS
  void dsymv_(const char* uplo, const int* n,
  	      const double* alpha, const double* a, const int* lda,
  	      const double* x, const int* incx, const double* beta,
  	      double* y, const int* incy);
  void dgemv_(const char* trans, const int* m, const int* n,
  	      const double* alpha, const double* a, const int* lda,
  	      const double* x, const int* incx, const double* beta,
  	      double* y, const int* incy);
  void daxpy_(const int *n, const double *da, const double *dx,
	      const int *incx, double *dy, const int *incy);
  double dnrm2_(const int *n, const double *dx, const int *incx);

  
} // end FORTRAN definitions

#endif
