//
// 
//                This source code is part of
// 
//                        M D S C T K
// 
//       Molecular Dynamics Spectral Clustering ToolKit
// 
//                        VERSION 1.2.5
// Written by Joshua L. Phillips.
// Copyright (c) 2012-2016, Joshua L. Phillips.
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
// official version at github.com/jlphillipsphd/mdsctk/.
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

// const double RAD2DEG = 180.0 / M_PI;
const size_t update_interval = 3;

CSC_matrix::CSC_matrix() {
  init();
}

CSC_matrix::CSC_matrix(const string filename) {
  init();
  ifstream csc;
  csc.open(filename.c_str());
  csc.read((char*) &n, (sizeof(int) / sizeof(char)));
  pcol = new int[n+1];
  csc.read((char*) pcol, (sizeof(int) / sizeof(char)) * (n+1));
  nnz = pcol[n];
  M = new double[nnz];
  irow = new int[nnz];
  csc.read((char*) irow, (sizeof(int) / sizeof(char)) * nnz);
  csc.read((char*) M, (sizeof(double) / sizeof(char)) * nnz);  
  if (csc.bad())
    cleanup();
  csc.close();
}

CSC_matrix::CSC_matrix(const CSC_matrix& Rhs) {
  init();
  copy(Rhs);
}

CSC_matrix::~CSC_matrix() {
  cleanup();
}

CSC_matrix& CSC_matrix::operator=(const CSC_matrix &Rhs) {
  if (this != &Rhs) {
    cleanup();
    copy(Rhs);
  }
  return *this;
}

void CSC_matrix::syslice(vector<int>& rc, CSC_matrix& csc) {
  vector<int> new_pcol;
  vector<int> new_irow;
  vector<double> new_M;
  
  int new_nnz = 0;
  for (int x = 0; x < rc.size(); x++) {
    new_pcol.push_back(new_nnz);
    // cout << "Starting new column (" << x << ") using " << rc[x] << endl;
    int z = 0;
    for (int y = pcol[rc[x]]; y < pcol[rc[x]+1]; y++) {
      for (; z < rc.size() && rc[z] <= irow[y]; z++) {
	// cout << "Working: "
	//      << x << "(" << rc[x] << ") "
	//      << y << "(" << irow[y] << ") "
	//      << z << "(" << rc[z] << ")";
	if (rc[z] == irow[y]) {
	  // cout << " [INS] " << z << " " << x;
	  new_irow.push_back(z);
	  new_M.push_back(M[y]);
	  new_nnz++;
	}
	// cout << endl;
      }
    }
  }
  new_pcol.push_back(new_nnz);
  
  // Construct the newbie...
  csc.cleanup();
  csc.n = rc.size();
  csc.nnz = new_nnz;
  csc.pcol = new int[new_pcol.size()];
  csc.irow = new int[new_irow.size()];
  csc.M = new double[new_M.size()];

  std::copy(new_pcol.begin(),new_pcol.end(),csc.pcol);
  std::copy(new_irow.begin(),new_irow.end(),csc.irow);
  std::copy(new_M.begin(),new_M.end(),csc.M);

  return;
}

void CSC_matrix::cleanup() {
  n=0;
  nnz=0;
  if (irow)
    delete [] irow;
  if (pcol)
    delete [] pcol;
  if (M)
    delete [] M;
  irow=NULL;
  pcol=NULL;
  M=NULL;
}

void CSC_matrix::init() {
  n=0;
  nnz=0;
  irow=NULL;
  pcol=NULL;
  M=NULL;
}

void CSC_matrix::copy(const CSC_matrix &Rhs) {
  n=Rhs.n;
  nnz=Rhs.nnz;
  pcol = new int[n+1];
  irow = new int[nnz];
  M = new double[nnz];
  memcpy(pcol,Rhs.pcol,sizeof(int)*(n+1));
  memcpy(irow,Rhs.irow,sizeof(int)*nnz);
  memcpy(M,Rhs.M,sizeof(double)*n);  
}

double& CSC_matrix::operator[](int x) {
  return M[x];
}

TOP_file::TOP_file(const string init_filename) : filename(init_filename),
						 natoms(0),
						 mass(NULL),
						 frame(NULL)
{
  read_topology();
  natoms = top.atoms.nr;
  if (natoms > 0) {
    mass = new ::real[natoms];
    for (int x = 0; x < natoms; x++)
      mass[x] = top.atoms.atom[x].m;
    center(frame);
  }
}

TOP_file::~TOP_file() {
  if (mass != NULL)
    delete [] mass;
}

int TOP_file::get_natoms() {
  return natoms;
}

::real* TOP_file::get_mass() {
  return mass;
}

coord_array TOP_file::get_frame_ptr() {
  return frame;
}

void TOP_file::center(coord_array frame) {
  reset_x(natoms,NULL,natoms,NULL,frame,mass);
}

::real TOP_file::rmsd(coord_array ref_frame, coord_array fit_frame) {
  do_fit(natoms,mass,ref_frame,fit_frame);
  return rmsdev(natoms,mass,ref_frame,fit_frame) * 10.0;
}

void TOP_file::com(coord_array frame, int n, int index[], rvec com) {
  gmx_calc_com(&top,frame,n,index,com);
}

void TOP_file::get_index(const string ndx_filename, int &ndx_n, int* &ndx_index, char* &ndx_group) {
  if (ndx_filename == "")
    ::get_index(&(top.atoms),NULL,1,&ndx_n,&ndx_index,&ndx_group);
  else
    ::get_index(&(top.atoms),ndx_filename.c_str(),1,&ndx_n,&ndx_index,&ndx_group);
}


void TOP_file::read_topology() {
  read_tps_conf(filename.c_str(), &top, &ePBC, &frame,
		NULL, box, TRUE);
}


XTC_file::XTC_file(const string init_filename) : filename(init_filename),
						 file(NULL),
						 step(1),
						 time(0.0),
						 prec(0.001),
						 bOK(1),
						 natoms(0),
						 frame(NULL)
{
  file = open_xtc(filename.c_str(),"r");
  read_first_xtc(file, &natoms, &step, &time, box, &frame, &prec, &bOK);
  close_xtc(file);
  if (natoms > 0) {
    file = open_xtc(filename.c_str(),"r");
    frame = new rvec[natoms];
  }
}

XTC_file::~XTC_file() {
  close_xtc(file); 
  delete [] frame;
}
  
int XTC_file::get_natoms() {
  return natoms;
}
float XTC_file::get_time() {
  return time;
}
float XTC_file::get_prec() {
  return prec;
}
int XTC_file::get_step() {
  return step;
}

coord_array XTC_file::get_next_frame_ptr() {
  if (natoms > 0 &&
      read_next_xtc(file, natoms, &step, &time, box, frame, &prec, &bOK))
    return frame;
  return NULL;
}

coord_array XTC_file::get_next_frame() {
  if (get_next_frame_ptr()) {
    coord_array old_frame = frame;
    frame = new rvec[natoms];
    return old_frame;
  }
  return NULL;
}

void copyright(const char* program_name) {

  cout << endl;
  cout << "   MDSCTK " << MDSCTK_VERSION_MAJOR << "." << MDSCTK_VERSION_MINOR;
  if (program_name)
    cout << " - " << program_name << endl;
  else
    cout << endl;
  cout << "   Copyright (C) 2013 Joshua L. Phillips" << endl;
  cout << "   MDSCTK comes with ABSOLUTELY NO WARRANTY; see LICENSE for details." << endl;
  cout << "   This is free software, and you are welcome to redistribute it" << endl;
  cout << "   under certain conditions; see README.md for details." << endl;
  cout << endl;

}

double getEPS() {
  double eps = 1.0;
  do { eps /= 2.0; } while (1.0 + (eps / 2.0) != 1.0);
  eps = sqrt(eps);
  return (eps);
}

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

double euclidean_distance(int size, double* reference, double* fitting) {
  double value = 0.0;
  for (int x = 0; x < size; x++)
    value += (reference[x] - fitting[x]) * (reference[x] - fitting[x]);
  return (sqrt(value));
}

double correlation_distance(int size, double* reference, double* fitting) {
  double rvalue = 0.0;
  double rsvalue = 0.0;
  double fvalue = 0.0;
  double fsvalue = 0.0;
  double dsize = (double) size;
  double value = 0.0;
  for (int x = 0; x < size; x++) {
    rvalue += reference[x];
    rsvalue += reference[x]*reference[x];
    fvalue += fitting[x];
    fsvalue += fitting[x]*fitting[x];
  }
  rsvalue = sqrt(((dsize*rsvalue)-(rvalue*rvalue))/(dsize*(dsize-1.0)));
  rvalue /= dsize;
  fsvalue = sqrt(((dsize*fsvalue)-(fvalue*fvalue))/(dsize*(dsize-1.0)));
  fvalue /= dsize;
  for (int x = 0; x < size; x++)
    value += (reference[x]-rvalue)*(fitting[x]-fvalue);
  value = (1.0 - (value / ((dsize - 1.0) * rsvalue * fsvalue))) / 2.0;
  if (value < 0.0)
    value = 0.0;
  return sqrt(value);
}

double euclidean_distance_sparse(int ref_size, int* ref_index, double* ref_data,
				 int fit_size, int* fit_index, double* fit_data) {
  double value = 0.0;
  int ref = 0;
  int fit = 0;

  for (ref = 0; ref < ref_size; ref++) {
    while (fit < fit_size && fit_index[fit] < ref_index[ref]) {
      value += (fit_data[fit] * fit_data[fit]);
      fit++;
    }
    if (fit < fit_size && ref_index[ref] == fit_index[fit]) {
      value += ((ref_data[ref] - fit_data[fit]) *
		(ref_data[ref] - fit_data[fit]));
      fit++;
    }
    else {
      value += (ref_data[ref] * ref_data[ref]);
    }
  }
  for (;fit < fit_size;fit++) {
    value += (fit_data[fit] * fit_data[fit]);
  }
  return (sqrt(value));
}

double entropic_affinity_sigma(vector<double> &A, int k, double b0,
			       double logK, double logN,
			       double B_lower, double B_upper) {
  int maxit = 20;
  double tol = 1e-10;
  double b = 0.0;
  double realmin = 2.225074e-308;

  // EPS
  double eps = 1.0;
  do { eps /= 2.0; } while (1.0 + (eps / 2.0) != 1.0);
  eps = sqrt(eps);

  if (b0<B_lower || b0>B_upper) {
    b = AVG(B_lower,B_upper);
  } else {
    b = b0;
  }

  int i = 1;
  vector<double> ed2;
  vector<double> m1v;
  ed2.resize(k);
  m1v.resize(k);
  double e = 0.0;
  double g = 0.0;
  double eg2 = 0.0;
  double m0 = 0.0;
  double m1 = 0.0;
  double m2 = 0.0;
  double m12 = 0.0;
  while (true) {
    double bE = exp(b);
    bool pbm = false;
    
    for (int x = 0; x < ed2.size(); x++)
      ed2[x] = exp(-SQR(A[x]) * bE);
    m0 = 0.0;
    for (int x = 0; x < ed2.size(); x++)
      m0 += ed2[x];

    if (m0<realmin) {
      e = -logK;
      pbm = true;
    } else {
      for (int x = 0; x < m1v.size(); x++)
	m1v[x] = ed2[x]*(SQR(A[x])/m0);
      m1 = 0.0;
      for (int x = 0; x < m1v.size(); x++)
	m1 += m1v[x];
      e = bE*m1 + log(m0) - logK;
    }

    // cout << i << " " << e << endl;

    if (fabs(e) < tol) { break; }

    // Very narrow bounds... no need to iterate
    if (B_upper-B_lower < 10.0*eps) { break; }

    if (e<0.0 && b<=B_upper) {
      B_upper = b;
    } else if (e>0.0 && b>=B_lower) {
      B_lower = b;
    }

    // Should check if e->Inf as well, but not implemented yet...
    pbm = pbm || e < -logK || e > logN-logK;
    
    if (!pbm) {
      if (i==maxit) {
	b = AVG(B_lower,B_upper);
	i=1;
	continue;
      }
      eg2 = SQR(bE);
      m2 = 0.0;
      for (int x = 0; x < m1v.size(); x++)
	m2 += m1v[x]*SQR(A[x]);
      m12 = SQR(m1)-m2;
      g = eg2*m12;
      if (g==0)
	pbm = true;
    }

    // Function problemd - bisect old bounds
    // Gradient problem - bisect new bounds
    if (pbm) {
      // Both problems? - just return
      double esqd_sum = 0.0;
      for (int x = 0; x < k; x++)
	esqd_sum += exp(-SQR(A[x]) * exp(B_lower)) + exp(-SQR(A[x]) * exp(B_upper));
      if (esqd_sum < 2.0*sqrt(realmin)) { break; }
      b = AVG(B_lower,B_upper);
      i = 1;
      continue;
    }

    // Newton step OK, update the bounds
    b += -e/g;
    // Out of bounds?
    if (b<B_lower || b > B_upper) {
      b = AVG(B_lower,B_upper);
      i = 0;
    }
    i++;
  }
  return (1.0 / sqrt(2.0 * exp(b)));
}

void entropic_affinity_sigmas(int n, int k, double K,
			      vector<double>* A, double* s) {
  double *B_lower = new double[n];
  double *B_upper = new double[n];
  int Ki = (int) ceil(K);
  double N = (double) k;
  double logK = log(K);
  double logN = log(N);
  double logNK = logN - logK;
  double p1 = 0.0;
  permutation<double> Kth;

  // EPS
  double eps = 1.0;
  do { eps /= 2.0; } while (1.0 + (eps / 2.0) != 1.0);
  eps = sqrt(eps);

  // Compute p1(N,logK)
  if (logK > log(sqrt(2.0*N))) {
    p1 = 3.0 / 4.0;
  }
  else {
    p1 = 1.0 / 4.0;
    for (int x = 0; x < 100; x++)
      p1 -= (-p1*log(p1/N)-logK)/(-log(p1/N)+1.0);
    p1 = 1.0 - (p1/2.0);
  }

  // NOTE there are no checks for ties. This
  // would be needed for stability...
  for (int x = 0; x < n; x++) {
    B_upper[x] = log((2.0*log(p1*(N-1.0)/(1.0-p1)))/
		(SQR(A[x][1])-SQR(A[x][0])));
    double bL1 = log((2.0*logNK/(1.0-(1.0/N)))/
		     (SQR(A[x][k-1])-SQR(A[x][0])));
    double bL2 = log((2.0*sqrt(logNK))/
		     sqrt((SQR(A[x][k-1])*SQR(A[x][k-1]))-
			  (SQR(A[x][0])*SQR(A[x][0]))));
    if (bL1 > bL2)
      B_lower[x] = bL1;
    else
      B_lower[x] = bL2;
  }


  // Sort
  Kth.data.resize(n);
  for (int x = 0; x < n; x ++) {
    Kth.data[x] = A[x][Ki-1];
  }
  Kth.sort();

  double b0 = AVG(B_lower[Kth.indices[0]],
		  B_upper[Kth.indices[0]]);
  // Setup loop
  for (vector<int>::iterator j = Kth.indices.begin(); j != Kth.indices.end(); j++) {
    s[*j] = entropic_affinity_sigma(A[*j],k,b0,
				    logK,logN,
				    B_lower[*j],B_upper[*j]);
    b0 = log((1.0/s[*j])*(1.0/s[*j])/2.0);
  }

  // for (int x = 0; x < n; x++)
  //   cout << x << " " << s[x] << endl;
    
  delete [] B_lower;
  delete [] B_upper;
}

int compare_edge(Db *db, const Dbt *key1, const Dbt *key2) {
  return compare_edge(db, key1, key2, NULL);
}

int compare_edge(Db *db, const Dbt *key1, const Dbt *key2, size_t* locp = NULL) {
  edge e1,e2;

  memcpy(&e1,key1->get_data(),sizeof(edge));
  memcpy(&e2,key2->get_data(),sizeof(edge));

  if (e1.from == e2.from)  {
    return (e1.to - e2.to);
  }
  return (e1.from - e2.from);
}

// Split edges
void split_edges(int current_index, Dbc *cursor, vector<int> &indices, vector<double> &distances) {

  edge myedge;
  double mydistance;
  Dbt key(&myedge,sizeof(edge));
  Dbt data(&mydistance,sizeof(double));
  key.set_ulen(sizeof(myedge));
  key.set_flags(DB_DBT_USERMEM);
  data.set_ulen(sizeof(double));
  data.set_flags(DB_DBT_USERMEM);

  indices.clear();
  distances.clear();
  
  if (cursor->get(&key, &data, DB_CURRENT) == 0) {
    do {
      if (myedge.from == current_index) {
	indices.push_back(myedge.to);
	distances.push_back(mydistance);
      }
      else
	break;
    } while (cursor->get(&key, &data, DB_NEXT) == 0);
  }
}

::real theta(::real pos1[], ::real pos2[], ::real pos3[],
	      bool degrees) {
  ::real L[3], Lnorm;
  ::real R[3], Rnorm;
  ::real angle;
  
  L[0] = pos1[0] - pos2[0];
  L[1] = pos1[1] - pos2[1];
  L[2] = pos1[2] - pos2[2];

  R[0] = pos3[0] - pos2[0];
  R[1] = pos3[1] - pos2[1];
  R[2] = pos3[2] - pos2[2];

  Lnorm = sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]);
  Rnorm = sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2]);
  angle = -(L[0]*R[0] + L[1]*R[1] + L[2]*R[2]) / (Lnorm * Rnorm);
  
  if (angle > 1.0) angle = 1.0;
  if (angle < -1.0) angle = -1.0;
  
  angle = acos( angle );
  if (degrees)
    angle = angle * RAD2DEG;
  
  return angle;
}

void crossprod(::real C[],
	       ::real x1, ::real y1, ::real z1,
	       ::real x2, ::real y2, ::real z2) {
  C[0] = ((y1 * z2) - (z1 * y2));
  C[1] = ((z1 * x2) - (x1 * z2));
  C[2] = ((x1 * y2) - (y1 * x2));
  return;
}

::real torsion(::real pos1[], ::real pos2[],
	       ::real pos3[], ::real pos4[], bool degrees) {
  ::real L[3], Lnorm;
  ::real R[3], Rnorm;
  ::real S[3];
  ::real angle;

  crossprod(L,
	    (pos2[0] - pos1[0]), (pos2[1] - pos1[1]), (pos2[2] - pos1[2]),
	    (pos3[0] - pos2[0]), (pos3[1] - pos2[1]), (pos3[2] - pos2[2]));
  crossprod(R,
	    (pos4[0] - pos3[0]), (pos4[1] - pos3[1]), (pos4[2] - pos3[2]),
	    (pos2[0] - pos3[0]), (pos2[1] - pos3[1]), (pos2[2] - pos3[2]));

  Lnorm = sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]);
  Rnorm = sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2]);

  crossprod(S, L[0], L[1], L[2], R[0], R[1], R[2]);
  angle = (L[0]*R[0] + L[1]*R[1] + L[2]*R[2]) / (Lnorm * Rnorm);

  if (angle > 1.0) angle = 1.0;
  if (angle < -1.0) angle = -1.0;

  angle = acos( angle );
  if (degrees)
    angle = angle * RAD2DEG;

  if ((S[0] * (pos3[0]-pos2[0]) + 
       S[1] * (pos3[1]-pos2[1]) +
       S[2] * (pos3[2]-pos2[2])) < 0 )
    angle = -angle;

  return angle;
}

void sample(int n, int k, int *sample) {
  int t = 0;
  int m = 0;
  double u = 0.0;
  while (m < k) {
    u = rand() / (RAND_MAX + 1.0);
    if ((n-t)*u >= k-m)
      t++;
    else
      sample[m++] = t++;	
  }
  return;
}

// Note FORTRAN convention...
void kmeans(int n, int d, int k, double *data, int *labels, int nstarts, int maxit) {
  int *kset = new int[k];
  double fit = (numeric_limits<double>::max)();
  int *clabels = new int[n];
  double *fits = new double[k];
  double *centers = new double[k*d];

  // cout << "Pointer " << kset << endl;

  // Preliminary initialization
  for (int x = 0; x < n; x++) {
    clabels[x] = 0;
  }
  for (int x = 0; x < k; x++) {
    fits[x] = (numeric_limits<double>::max)();
    for (int y = 0; y < d; y++) {
      centers[(x*d)+y] = 0.0;
    }
  }

  // cout << "Preliminary startup..." << endl;

  // Start computation
  for (int start = 0; start < nstarts; start++) {
    
    // Initialize run
    // cout << "Starting run #" << start << endl;

    // Get random centers...
    sample(n,k,kset);
    for (int x = 0; x < k; x++) {
      for (int y = 0; y < d; y++) 
	centers[(x*d)+y] = data[(y*n)+kset[x]];
      fits[x] = 0.0;
    }

    // cout << "Initial centers: " << endl;
    // for (int x = 0; x < k; x++) {
    //   cout << x << " (" << kset[x] << ") ";
    //   for (int y = 0; y < d; y++)
    // 	cout << centers[(x*d)+y] << " ";
    //   cout << endl;
    // }

    // Calc fits
    for (int x = 0; x < n; x++) {
      int currenti = -1;
      double currentd = (numeric_limits<double>::max)();
      for (int y = 0; y < k; y++) {
	double dist = 0.0;
	for (int z = 0; z < d; z++) {
	  dist += SQR((data[(z*n)+x] - centers[(y*d)+z]));
	}
	if (dist < currentd) {
	  currenti = y;
	  currentd = dist;
	}
      }
      // cout << "Pt: " << x << "( ";
      // for (int y = 0; y < d; y++)
      // 	cout << data[(y*n)+x] << " ";
      // cout << ") " << currenti << " " << currentd << endl;
      fits[currenti] += currentd;
      clabels[x] = currenti;
    }

    // cout << "Initial assignment: " << endl;
    // for (int x = 0; x < n; x++)
    //   cout << clabels[x] << " ";
    // cout << endl;

    // Iterate
    for (int it = 0; it < maxit; it++) {

      // Calc new centers
      for (int x = 0; x < k; x++) {
	kset[x] = 0;
	fits[x] = 0.0;
	for (int y = 0; y < d; y++)
	  centers[(x*d)+y] = 0.0;
      }
      for (int x = 0; x < n; x++) {
	for (int y = 0; y < d; y++) {
	  centers[(clabels[x]*d)+y] += data[(y*n)+x];
	}
	kset[clabels[x]]++;
      }
      for (int x = 0; x < k; x++)
	for (int y = 0; y < d; y++)
	  centers[(x*d)+y] /= (double) kset[x];

      // cout << "New centers: " << endl;
      // for (int x = 0; x < k; x++) {
      // 	cout << x << " ";
      // 	for (int y = 0; y < d; y++)
      // 	  cout << centers[(x*d)+y] << " ";
      // 	cout << endl;
      // }
      
      // Compute change in assignment -- if any...
      bool change = false;
      for (int x = 0; x < n; x++) {
	int currenti = -1;
	double currentd = (numeric_limits<double>::max)();
	for (int y = 0; y < k; y++) {
	  double dist = 0.0;
	  for (int z = 0; z < d; z++) {
	    dist += SQR((data[(z*n)+x] - centers[(y*d)+z]));
	  }
	  if (dist < currentd) {
	    currenti = y;
	    currentd = dist;
	  }
	}
	if (currenti != clabels[x])
	  change = true;
	fits[currenti] += currentd;
	clabels[x] = currenti;
      }

      // double ss = 0.0;
      // for (int x = 0; x < k; x++)
      // 	ss += fits[x];
      // cout << "Iteration " << it << ": " << ss << endl;

      if (!change) {
	// cout << "Terminating on iteration " << it << endl;
	break;
      }
      
    } // end it

    // cout << "Final assignment: " << endl;
    // for (int x = 0; x < n; x++)
    //   cout << clabels[x] << " ";
    // cout << endl;

    // Update current partition if better than others...
    double ss = 0.0;
    for (int x = 0; x < k; x++)
      ss += fits[x];
    if (ss < fit) {
      fit = ss;
      memcpy(labels,clabels,sizeof(int)*n);
    }
  }

  // cout << "Pointer " << kset << endl;

  delete [] kset;
  delete [] clabels;
  delete [] fits;
  delete [] centers;
  return;
}

vector<int> select(vector<int> &set, int k, int *labels) {
  vector<int> myset;
  for (int x = 0; x < set.size(); x++)
    if (labels[x] == k)
      myset.push_back(set[x]);
  return myset;
}

int runARPACK(int nev, CSC_matrix &A, double* &d, double* &Z) {
  // ARPACK variables...
  int ido = 0;
  char bmat = 'I';
  char which[2];
  which[0] = 'L';
  which[1] = 'A';
  double tol = 0.0;
  double *resid = new double[A.n];
  // NOTE: Need about one order of magnitude more arnoldi vectors to
  // converge for the normalized Laplacian (according to residuals...)
  int ncv = ((10*nev+1)>A.n)?A.n:(10*nev+1);
  double *V = new double[(ncv*A.n)+1];
  int ldv = A.n;
  int *iparam = new int[12];
  iparam[1] = 1;
  iparam[3] = 100 * nev;
  iparam[4] = 1;
  iparam[7] = 1;
  int *ipntr = new int[15];
  double *workd = new double[(3*A.n)+1];
  int lworkl = ncv*(ncv+9);
  double *workl = new double[lworkl+1];
  int info = 0;
  int rvec = 1;
  char HowMny = 'A';
  int *lselect = new int[ncv];
  d = new double[nev];
  // Z = &V[1];
  Z = V;
  int ldz = A.n;
  double sigma = 0.0;
 
  while (ido != 99) {
    dsaupd_(&ido, &bmat, &A.n, which,
	    &nev, &tol, resid,
	    &ncv, &V[1], &ldv,
	    &iparam[1], &ipntr[1], &workd[1],
	    &workl[1], &lworkl, &info);
    
    if (ido == -1 || ido == 1) {
      // Matrix-vector multiplication
      sp_dsymv(A.n,A.irow,A.pcol,A.M,
	       &workd[ipntr[1]],
	       &workd[ipntr[2]]);
    }
  }
    
  dseupd_(&rvec, &HowMny, lselect,
	  d, Z, &ldz,
	  &sigma, &bmat, &A.n,
	  which, &nev, &tol,
	  resid, &ncv, &V[1],
	  &ldv, &iparam[1], &ipntr[1],
	  &workd[1], &workl[1],
	  &lworkl, &info);

  ido = iparam[5];

  delete [] resid;
  delete [] lselect;
  delete [] iparam;
  delete [] ipntr;
  delete [] workd;
  delete [] workl;

  return ido;
}

int runARPACK2(int nev, CSC_matrix &A, double* &d, double* &Z) {
  // ARPACK variables...
  int ido = 0;
  char bmat = 'I';
  char which[2];
  which[0] = 'L';
  which[1] = 'R';
  double tol = 0.0;
  double *resid = new double[A.n];
  // NOTE: Need about one order of magnitude more arnoldi vectors to
  // converge for the normalized Laplacian (according to residuals...)
  int ncv = ((10*nev+1)>A.n)?A.n:(10*nev+1);
  double *V = new double[(ncv*A.n)+1];
  int ldv = A.n;
  int *iparam = new int[12];
  iparam[1] = 1;
  iparam[3] = 100 * nev;
  iparam[4] = 1;
  iparam[7] = 1;
  int *ipntr = new int[15];
  double *workd = new double[(3*A.n)+1];
  int lworkl = ncv*(3*ncv+9);
  double *workl = new double[lworkl+1];
  int info = 0;
  int rvec = 1;
  char HowMny = 'A';
  int *lselect = new int[ncv];
  d = new double[nev+1];
  double *di = new double[nev+1];
  // Z = &V[1];
  Z = V;
  int ldz = A.n;
  double sigmar = 0.0;
  double sigmai = 0.0;
  double *workev = new double[(3*ncv)+1];

  while (ido != 99) {
    dnaupd_(&ido, &bmat, &A.n, which,
	    &nev, &tol, resid,
	    &ncv, &V[1], &ldv,
	    &iparam[1], &ipntr[1], &workd[1],
	    &workl[1], &lworkl, &info);
    
    if (ido == -1 || ido == 1) {
      // Matrix-vector multiplication
      sp_dgemv(A.n,A.irow,A.pcol,A.M,
	       &workd[ipntr[1]],
	       &workd[ipntr[2]]);
    }
  }
    
  dneupd_(&rvec, &HowMny, lselect,
	  d, di, Z, &ldz,
	  &sigmar, &sigmai, &workev[1],
	  &bmat, &A.n,
	  which, &nev, &tol,
	  resid, &ncv, &V[1],
	  &ldv, &iparam[1], &ipntr[1],
	  &workd[1], &workl[1],
	  &lworkl, &info);

  ido = iparam[5];

  cout << "Number of eigenvalues/vectors found: " << ido << endl;

  cout << "Resid:" << endl;
  for (int x = 0; x < nev; x++)
    cout << resid[x] << " ";
  cout << endl;

  for (int x = nev-1; x >= 0; x--) {
    cout << "Eigval: " << d[x] << " +i" << di[x] << " (" << (1.0/(1.0-d[x])) << ")" << endl;
  }

  ofstream out;
  out.open("eigenvectors.dat");
  for (int x = nev-1; x >= 0; x--) {
    for (int y = 0; y < 2*A.n; y++)
      out << Z[2*(x*A.n)+y] << " ";
  out << endl;
  }
  out.close();

  delete [] di;
  delete [] resid;
  delete [] lselect;
  delete [] iparam;
  delete [] ipntr;
  delete [] workd;
  delete [] workl;
  delete [] workev;

  return ido;
}

int gmx_calc_com(t_topology *top, rvec x[], int nrefat, int index[], rvec xout)
{
  float                mass, total_mass;
  
  if (!top) {
    cerr << "No masses available but mass weighting was requested [COM CALCULATION FAILURE!]" << endl;
    return -1;
  }
  xout[0] = xout[1] = xout[2] = 0.0;
  total_mass = 0.0;
  for (int i = 0; i < nrefat; i++) {
    int ai = index[i];
    mass = top->atoms.atom[ai].m;
    for (int j = 0; j < 3; ++j) {
      xout[j] += mass * x[ai][j];
    }
    total_mass += mass;
  }
  xout[0] /= total_mass;
  xout[1] /= total_mass;
  xout[2] /= total_mass;
  return 0;
}


#if (GMX_VERSION < 20160000)
gmx_bool read_tps_conf(const char *infile, t_topology *top, int *ePBC,
                       rvec **x, rvec **v, matrix box, gmx_bool requireMasses) {
  char buf[256];
  return read_tps_conf(infile, buf, top, ePBC, x, v, box, requireMasses);
}
#endif
