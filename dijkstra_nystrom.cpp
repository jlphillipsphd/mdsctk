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

// C
#include <stdio.h>

// C++
#include <fstream>
#include <iostream>

// BOOST C++
#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

// OpenMP
#include <omp.h>

// Local
#include "config.h"

#define INDEX(SIZE,I,J) ((SIZE * J) + I)

using namespace boost;

int main(int argc, char* argv[]) {

  if (argc != 2) {
    std::cerr << std::endl;
    std::cerr << "   MDSCTK " << MDSCTK_VERSION_MAJOR << "." << MDSCTK_VERSION_MINOR << std::endl;
    std::cerr << std::endl;
    std::cerr << "Usage: " << argv[0] << " [#threads]" << std::endl;
    std::cerr << "   Computes the shortest paths through the provided CSC sparse matrix." << std::endl;
    std::cerr << "   NOTE: Implicit zeros are assumed to be edges of INF weight, but all" << std::endl;
    std::cerr << "   diagonal entries are assumed zero (even if they are set to non-zero!)." << std::endl;
    std::cerr << std::endl;
    return -1;
  }

  // General
  int     n;   // Dimension of the problem.
  int     m;   // Outer dimension
  int     nthreads;  // # threads

  // Main affinity matrix
  int     nnzA;
  int     *irowA;
  int     *pcolA;
  int     nbA;  // Block size for reading data
  double  *A;   // Pointer to an array that stores the lower
		// triangular elements of A.

  // Expanded affinity matrix
  int     nnzB;
  int     *irowB;
  int     *pcolB;
  int     nbB;  // Block size for reading data
  double  *B;   // Pointer to an array that stores the
		// sparse elements of B.

  // Boost data types
  typedef adjacency_list < listS, vecS, directedS,
    no_property, property < edge_weight_t, double > > graph_t;
  typedef graph_traits < graph_t >::vertex_descriptor vertex_descriptor;
  typedef graph_traits < graph_t >::edge_descriptor edge_descriptor;
  typedef std::pair<int, int> Edge;

  // File input streams
  std::ifstream distances;
  std::ifstream row_pointers;
  std::ifstream col_pointers;

  // File output stream
  std::ofstream apsp;
  std::ofstream alsp;

  // Open output files
  apsp.open("apsp.dat");
  alsp.open("alsp.dat");
  nthreads=atoi(argv[1]);
  omp_set_num_threads(nthreads);

  // Read sparse matrix data
  distances.open("sym_distances.dat");
  col_pointers.open("sym_col_indices.dat");
  row_pointers.open("sym_row_indices.dat");
 
  if (col_pointers.bad()) {
    std::cerr << "*** ERROR ***" << std::endl;
    std::cerr << "   Could not open column pointer file: " << argv[1] << std::endl;
    return -1;
  }

  if (row_pointers.bad()) {
    std::cerr << "*** ERROR ***" << std::endl;
    std::cerr << "   Could not open row pointer file: " << argv[2] << std::endl;
    return -1;
  }

  if (distances.bad()) {
    std::cerr << "*** ERROR ***" << std::endl;
    std::cerr << "   Could not open symmetric CSC file: " << argv[3] << std::endl;
    return -1;
  }

  std::cerr << "Reading landmark sparse matrix data...";

  // Determine size of column pointers vector
  col_pointers.seekg(0,std::ios::end);
  n = (col_pointers.tellg() * sizeof(char) / sizeof(int));
  col_pointers.seekg(0,std::ios::beg);
  pcolA = new int[n];
  col_pointers.read((char*) pcolA, (sizeof(int) / sizeof(char)) * n);
  n--;
  nnzA = pcolA[n];
  nbA = 0;
  for (int x = 1; x <= n; x++)
    if ((pcolA[x] - pcolA[x-1]) > nbA)
      nbA = pcolA[x] - pcolA[x-1];
  A = new double[nbA];
  irowA = new int[nbA];

  // Read row index and matrix entries (a)
  graph_t g(n);
  property_map<graph_t, edge_weight_t>::type weightmap = get(edge_weight, g);
  for (int i = 0; i < n; i++) {
    nbA = pcolA[i+1] - pcolA[i];
    row_pointers.read((char*) irowA, (sizeof(int) / sizeof(char)) * nbA);
    distances.read((char*) A, (sizeof(double) / sizeof(char)) * nbA);
    for (int j = 0; j < nbA; j++) {
      edge_descriptor e; bool inserted;
      tie(e, inserted) = add_edge(i, irowA[j], g);
      weightmap[e] = A[j];
      tie(e, inserted) = add_edge(irowA[j], i, g);
      weightmap[e] = A[j];
    }
  }

  std::cerr << "done." <<std::endl;

  distances.close();
  col_pointers.close();
  row_pointers.close();

  std::vector<vertex_descriptor> *p[n];
  std::vector<double> *d[n];
  for (int j = 0; j < n; j++) {
    p[j] = new std::vector<vertex_descriptor>(num_vertices(g));
    d[j] = new std::vector<double>(num_vertices(g));
  }

  std::cerr << "Total number of vertices: " << num_vertices(g) << std::endl;
  std::cerr << "Total number of edges (symmetric): " << num_edges(g) / 2 << std::endl;
  std::cerr << "Computing shortest paths using Dijkstra's algorithm..." << std::endl;
  std::cerr << "Complete: 0%";

#pragma omp parallel for
  for (int j = 0; j < n; j++) {
    
    if (omp_get_thread_num() == 0)
      std::cerr << "\rComplete: " << (100 * j * nthreads) / n << "%";
    
    vertex_descriptor s = vertex(j, g);
    dijkstra_shortest_paths(g, s, &(p[j]->at(0)), &(d[j]->at(0)), weightmap, get(vertex_index, g),
			    std::less<double>(), closed_plus<double>(),
			    (std::numeric_limits<double>::max)(), 0.0,
			    default_dijkstra_visitor());
    
  }
  
  std::cerr << "\rComplete: 100%" << std::endl;
  
  for (int j = 0; j < n; j++) {
    apsp.write((char*) &(d[j]->at(j)), (sizeof(double) / sizeof(char)) * (n-j));
    delete p[j];
    delete d[j];
  }

  delete [] pcolA;
  delete [] irowA;
  delete [] A;

  // delete [] pcolB;
  // delete [] irowB;
  // delete [] B;

  return 0;
}
