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

int main(int argc, char* argv[]) {

  const char* program_name = "dijkstra";
  bool optsOK = true;
  gmx::initForCommandLine(&argc,&argv);
  copyright(program_name);
  cout << "   Computes the shortest paths through the provided CSC sparse matrix." << endl;
  cout << "   NOTE: Implicit zeros are assumed to be edges of INF weight, but all" << endl;
  cout << "   diagonal entries are assumed zero (even if they are set to non-zero!)." << endl;
  cout << endl;
  cout << "   Use -h or --help to see the complete list of options." << endl;
  cout << endl;

  // Option vars...
  int nthreads = 0;
  string ssm_filename;
  string o_filename;

  // Declare the supported options.
  po::options_description cmdline_options;
  po::options_description program_options("Program options");
  program_options.add_options()
    ("help,h", "show this help message and exit")
    ("threads,t", po::value<int>(&nthreads)->default_value(omp_get_max_threads()>omp_get_num_procs()?omp_get_num_procs():omp_get_max_threads()), "Input:  Number of threads to start (int)")
    ("ssm-file,f", po::value<string>(&ssm_filename)->default_value("distances.ssm"), "Input:  Symmetric CSC matrix file (string:filename)")
    ("output-file,o", po::value<string>(&o_filename)->default_value("apsp.dat"), "Output: All pairs shortest paths (string:filename)")
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
  cout << "threads =     " << nthreads << endl;
  cout << "ssm-file =    " << ssm_filename << endl;
  cout << "output-file = " << o_filename << endl;
  cout << endl;

  // Boost data types
  typedef boost::adjacency_list < boost::listS, boost::vecS, boost::directedS,
				  boost::no_property, boost::property < boost::edge_weight_t, double > > graph_t;
  typedef boost::graph_traits < graph_t >::vertex_descriptor vertex_descriptor;
  typedef boost::graph_traits < graph_t >::edge_descriptor edge_descriptor;
  typedef pair<int, int> Edge;

  // File output stream
  ofstream apsp;

  // Read sparse matrix data
  apsp.open(o_filename.c_str());
  omp_set_num_threads(nthreads);

  cout << "Reading sparse matrix data...";
  CSC_matrix *A = new CSC_matrix(ssm_filename.c_str());
  int n = A->n;

  // Read row index and matrix entries (a)
  graph_t g(n);
  boost::property_map<graph_t, boost::edge_weight_t>::type weightmap = get(boost::edge_weight, g);
  for (int i = 0; i < n; i++) {
    for (int j = A->pcol[i]; j < A->pcol[i+1]; j++) {
      edge_descriptor e; bool inserted;
      tie(e, inserted) = add_edge(i, A->irow[j], g);
      weightmap[e] = A->M[j];
      tie(e, inserted) = add_edge(A->irow[j], i, g);
      weightmap[e] = A->M[j];
    }
  }
  delete A;
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

  // Timer for ETA
  time_t start = std::time(0);
  time_t last = start;

#pragma omp parallel for
  for (int j = 0; j < n; j++) {
    
    // Update user of progress
    if (std::time(0) - last > update_interval) {
      last = std::time(0);
      time_t eta = start + ((last-start) * n / j);
      cout << "\rFrame: " << j << ", will finish " 
	   << string(std::ctime(&eta)).substr(0,20);
      cout.flush();
    }
    
    vertex_descriptor s = vertex(j, g);
    boost::dijkstra_shortest_paths(g, s, &(p[j]->at(0)),
				   &(d[j]->at(0)), weightmap,
				   get(boost::vertex_index, g),
				   less<double>(), boost::closed_plus<double>(),
				   (numeric_limits<double>::max)(), 0.0,
				   boost::default_dijkstra_visitor());
    
  }
  
  cout << endl << endl;
  
  for (int j = 0; j < n; j++) {
    apsp.write((char*) &(d[j]->at(j)), (sizeof(double) / sizeof(char)) * (n-j));
    delete p[j];
    delete d[j];
  }

  return 0;
}
