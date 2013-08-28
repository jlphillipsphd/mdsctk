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

// OpenMP
#include <omp.h>

// Local
#include "config.h"
#include "mdsctk.h"

int main(int argc, char* argv[]) {

  if (argc != 2) {
    cout << endl;
    cout << "   MDSCTK " << MDSCTK_VERSION_MAJOR << "." << MDSCTK_VERSION_MINOR << endl;
    cout << endl;
    cout << "Usage: " << argv[0] << " [#threads]" << endl;
    cout << "   Computes the shortest paths through the provided CSC sparse matrix." << endl;
    cout << "   NOTE: Implicit zeros are assumed to be edges of INF weight, but all" << endl;
    cout << "   diagonal entries are assumed zero (even if they are set to non-zero!)." << endl;
    cout << endl;
    return -1;
  }

  int     n;   // Dimension of the problem.
  int     nnz; // Number of non-zero elements
  int     *irow; // CSC Row indices (nnz)
  int     *pcol; // CSC Col indices (n+1)
  double  *a;   // CSC LT Matrix
  int     nb;  // Block size for reading data
  int     nthreads;  // # threads

  // Boost data types
  typedef boost::adjacency_list < listS, vecS, directedS,
			   no_property, property < edge_weight_t, double > > graph_t;
  typedef boost::graph_traits < graph_t >::vertex_descriptor vertex_descriptor;
  typedef boost::graph_traits < graph_t >::edge_descriptor edge_descriptor;
  typedef pair<int, int> Edge;

  // File input streams
  ifstream distances;
  ifstream row_pointers;
  ifstream col_pointers;

  // File output stream
  ofstream apsp;

  // Read sparse matrix data
  distances.open("sym_distances.dat");
  col_pointers.open("sym_col_indices.dat");
  row_pointers.open("sym_row_indices.dat");
  apsp.open("apsp.dat");
  nthreads=atoi(argv[1]);
  omp_set_num_threads(nthreads);

  if (col_pointers.bad()) {
    cout << "*** ERROR ***" << endl;
    cout << "   Could not open column pointer file: " << argv[1] << endl;
    return -1;
  }

  if (row_pointers.bad()) {
    cout << "*** ERROR ***" << endl;
    cout << "   Could not open row pointer file: " << argv[2] << endl;
    return -1;
  }

  if (distances.bad()) {
    cout << "*** ERROR ***" << endl;
    cout << "   Could not open symmetric CSC file: " << argv[3] << endl;
    return -1;
  }

  cout << "Reading sparse matrix data...";

  // Determine size of column pointers vector
  col_pointers.seekg(0,ios::end);
  n = (col_pointers.tellg() * sizeof(char) / sizeof(int));

  // Read column pointers
  col_pointers.seekg(0,ios::beg);
  pcol = new int[n];
  col_pointers.read((char*) pcol, (sizeof(int) / sizeof(char)) * n);
  n--;
  nnz = pcol[n];
  nb = 0;
  for (int x = 1; x <= n; x++)
    if ((pcol[x] - pcol[x-1]) > nb)
      nb = pcol[x] - pcol[x-1];
  a = new double[nb];
  irow = new int[nb];

  // Read row index and matrix entries (a)
  graph_t g(n);
  boost::property_map<graph_t, edge_weight_t>::type weightmap = get(edge_weight, g);
  for (int i = 0; i < n; i++) {
    nb = pcol[i+1] - pcol[i];
    row_pointers.read((char*) irow, (sizeof(int) / sizeof(char)) * nb);
    distances.read((char*) a, (sizeof(double) / sizeof(char)) * nb);
    for (int j = 0; j < nb; j++) {
      edge_descriptor e; bool inserted;
      tie(e, inserted) = add_edge(i, irow[j], g);
      weightmap[e] = a[j];
      tie(e, inserted) = add_edge(irow[j], i, g);
      weightmap[e] = a[j];
    }
  }

  cout << "done." <<endl;

  vector<vertex_descriptor> *p[n];
  vector<double> *d[n];
  for (int j = 0; j < n; j++) {
    p[j] = new vector<vertex_descriptor>(num_vertices(g));
    d[j] = new vector<double>(num_vertices(g));
  }

  cout << "Total number of vertices: " << num_vertices(g) << endl;
  cout << "Total number of edges (symmetric): " << num_edges(g) / 2 << endl;
  cout << "Computing shortest paths using Dijkstra's algorithm..." << endl;
  cout << "Complete: 0%";

#pragma omp parallel for
  for (int j = 0; j < n; j++) {
    
    if (omp_get_thread_num() == 0)
      cout << "\rComplete: " << (100 * j * nthreads) / n << "%";
    
    vertex_descriptor s = vertex(j, g);
    boost::dijkstra_shortest_paths(g, s, &(p[j]->at(0)),
				   &(d[j]->at(0)), weightmap,
				   get(vertex_index, g),
				   less<double>(), closed_plus<double>(),
				   (numeric_limits<double>::max)(), 0.0,
				   boost::default_dijkstra_visitor());
    
  }
  
  cout << "\rComplete: 100%" << endl;
  
  for (int j = 0; j < n; j++) {
    apsp.write((char*) &(d[j]->at(j)), (sizeof(double) / sizeof(char)) * (n-j));
    delete p[j];
    delete d[j];
  }

  delete [] pcol;
  delete [] irow;
  delete [] a;

  return 0;
}
