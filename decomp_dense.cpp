// Local
#include "config.h"
#include "mdsctk.h"

/* Main program */
int main() {
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
  distances.open("apsp.dat");
  eigenvalues.open("eigenvalues.dat");
  eigenvectors.open("eigenvectors.dat");

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
