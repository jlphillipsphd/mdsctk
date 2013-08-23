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

// Standard
// C
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
// C++
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <algorithm>

// Boost
#include <boost/program_options.hpp>

// Berkeley DB
#include <db_cxx.h>

// Local
#include "config.h"
#include "mdsctk.h"

namespace po = boost::program_options;
using namespace std;

// Edges struct, comparison, and sorting routines
struct edge {
  int from;
  int to;  
};

int compare_edge(Db *db, const Dbt *key1, const Dbt *key2) {
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

int main(int argc, char *argv[]) {

  const char* program_name = "make_gesparse";
  bool optsOK = true;
  copyright(program_name);
  cout << "   Converts the results from knn_* into general CSC format." << endl;
  cout << endl;
  cout << "   Normally, the number of nearest neighbors in the input" << endl;
  cout << "   distances is used for constructing the CSC matrix." << endl;
  cout << "   However, you can set output-knn <= knn in order to" << endl;
  cout << "   subselect the number of neighbors to consider in the" << endl;
  cout << "   CSC representation. This makes it easy to store a" << endl;
  cout << "   large number of neighbors using knn_* but then use" << endl;
  cout << "   a subset for, say, computing approximate geodesic" << endl;
  cout << "   distances." << endl;
  cout << endl;
  cout << "   Use -h or --help to see the complete list of options." << endl;
  cout << endl;

  // Option vars...
  int k;
  int maxk;
  string i_filename;
  string d_filename;
  string o_filename;

  // Declare the supported options.
  po::options_description cmdline_options;
  po::options_description program_options("Program options");
  program_options.add_options()
    ("help,h", "show this help message and exit")
    ("knn,k", po::value<int>(&maxk), "Input:  K-nearest neighbors (int)")
    ("output-knn,n", po::value<int>(&k), "Input:  K-nn to keep in output (int)")
    ("index-file,i", po::value<string>(&i_filename)->default_value("indices.dat"), "Input:  Index file (string:filename)")
    ("distance-file,d", po::value<string>(&d_filename)->default_value("distances.dat"), "Input:  Distances file (string:filename)")
    ("output-file,o", po::value<string>(&o_filename)->default_value("distances.gsm"), "Output: General sparse matrix file (string:filename)")    
    ;
  cmdline_options.add(program_options);

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
  po::notify(vm);    

  if (vm.count("help")) {
    cout << "usage: " << program_name << " [options]" << endl;
    cout << endl;
    cout << cmdline_options << endl;
    return 1;
  }
  if (!vm.count("knn")) {
    cout << "ERROR: --knn not supplied." << endl;
    cout << endl;
    optsOK = false;
  }
  if (!vm.count("output-knn"))
    k = maxk;

  // Check
  if (maxk < k) {
    cout << "ERROR: Output k (" << k << ") is not less than the input k (" << maxk << ")." << endl;
    cout << endl;
    optsOK = false;
  }

  if (!optsOK) {
    return -1;
  }

  cout << "Running with the following options:" << endl;
  cout << "knn =            " << maxk << endl;
  cout << "output-knn =     " << k << endl;
  cout << "index-file =     " << i_filename << endl;
  cout << "distances-file = " << d_filename << endl;
  cout << "output-file =    " << o_filename << endl;
  cout << endl;

  int current_index = 0;
  int *current_indices;
  double *current_distances;
  vector<int> *sorted_indices;
  vector<double> *sorted_distances;

  stringstream nsd_filename;
  nsd_filename << "." << getpid() << ".nsd";
  stringstream nsr_filename;
  nsr_filename << "." << getpid() << ".nsr";
  stringstream nsc_filename;
  nsc_filename << "." << getpid() << ".nsc";

  int nonsym_col = 0;
  current_indices = new int[maxk];
  current_distances = new double[maxk];
  sorted_indices = new vector<int>;
  sorted_distances = new vector<double>;

  ifstream distances;
  ifstream indices;
  ofstream nonsym_distances;
  ofstream nonsym_row_indices;
  ofstream nonsym_col_indices;
  ofstream o_file;

  // Open files for reading/writing
  distances.open(d_filename.c_str(),ios::in | ios::binary);
  if (!distances.good()) {
    cout << "***ERROR***" << endl;
    cout << "Could not open file: " << d_filename << endl;
    cout << endl;
    return -1;
  }
  indices.open(i_filename.c_str(),ios::in | ios::binary);
  if (!indices.good()) {
    cout << "***ERROR***" << endl;
    cout << "Could not open file: " << i_filename << endl;
    cout << endl;
    return -1;
  }
  nonsym_distances.open(nsd_filename.str().c_str(),ios::out | ios::binary);
  if (!nonsym_distances.good()) {
    cout << "***ERROR***" << endl;
    cout << "Could not open file: " << nsd_filename << endl;
    cout << endl;
    return -1;
  }
  nonsym_row_indices.open(nsr_filename.str().c_str(),ios::out | ios::binary);
  if (!nonsym_row_indices.good()) {
    cout << "***ERROR***" << endl;
    cout << "Could not open file: " << nsr_filename << endl;
    cout << endl;
    return -1;
  }
  nonsym_col_indices.open(nsc_filename.str().c_str(),ios::out | ios::binary);
  if (!nonsym_col_indices.good()) {
    cout << "***ERROR***" << endl;
    cout << "Could not open file: " << nsc_filename << endl;
    cout << endl;
    return -1;
  }
  o_file.open(o_filename.c_str(),ios::out | ios::binary);
  if (!o_file.good()) {
    cout << "***ERROR***" << endl;
    cout << "Could not open file: " << o_filename << endl;
    cout << endl;
    return -1;
  }

  // Setup the database
  Db db(NULL, 0);               // Instantiate the Db object
  Dbc *cursor = NULL;
  char db_file_name[100];
  edge myedge;
  double mydistance;
  int index;
      
  Dbt key(&myedge,sizeof(edge));
  Dbt data(&mydistance,sizeof(double));
  key.set_ulen(sizeof(myedge));
  key.set_flags(DB_DBT_USERMEM);
  data.set_ulen(sizeof(double));
  data.set_flags(DB_DBT_USERMEM);
  
  sprintf(db_file_name,"nonsym_%d.db",getpid());
  // Initialize database
  // Set up the btree comparison function for this database
  db.set_bt_compare(compare_edge);
  u_int32_t oFlags = DB_CREATE; // Open flags;
  // | this with DB_RDONLY to make the connection read-only
  // Set to 0 to just open an existing database...
  // Open the database
  db.open(NULL,                // Transaction pointer 
	  db_file_name,          // Database file name 
	  NULL,                // Optional logical database name
	  DB_BTREE,            // Database access method
	  oFlags,              // Open flags
	  0);                  // File mode (using defaults)

  // Read a set of distances and indices
  distances.read((char*) current_distances, sizeof(double) * maxk);
  indices.read((char*) current_indices, sizeof(int) * maxk);

  while (!distances.eof() || !indices.eof()) {
 
    if (current_index % 1000 == 0)
      cout << "Frame number: " << current_index << "\r";
    // cout << current_index << " : ";
    
    // Print basic information
    // cout << "Starting work on " << current_index << ":";
    // for (int x = 0; x < k; x++)
    //   cout << " " << current_indices[x];
    // cout << endl;

    for (int x = 0; x < k; x++) {
            
      // Do Work
      myedge.from = current_index;	
      myedge.to = current_indices[x];
      mydistance = current_distances[x];
      if (db.put(NULL, &key, &data, 0) != 0) {
	cout << "Could not insert edge: "
	     << myedge.from << " "
	     << myedge.to << " -> "
	     << mydistance << endl;
      }

      // This is a non-symmetric system, so no double entries...
      // myedge.to = current_index;
      // myedge.from = current_indices[x];
      // mydistance = current_distances[x];
      // if (db.put(NULL, &key, &data, 0) != 0) {
      // 	cout << "Ccould not insert edge: "
      // 	     << myedge.from << " "
      // 	     << myedge.to << " -> "
      // 	     << mydistance << endl;
      // }
    }

    // Read a set of distances and indices
    distances.read((char*) current_distances, sizeof(double) * maxk);
    indices.read((char*) current_indices, sizeof(int) * maxk);
    current_index++;
  }

  // Initialize cursor
  db.cursor(NULL, &cursor, 0); 
  cursor->get(&key, &data, DB_FIRST);
    
  for (int col = 0; col < current_index; col++) {
    
    split_edges(col, cursor, *sorted_indices, *sorted_distances);
    
    nonsym_col_indices.write((char*) &nonsym_col, sizeof(int));
    nonsym_row_indices.write((char*) &(sorted_indices->front()), 
			  sizeof(int) * sorted_indices->size());
    nonsym_distances.write((char*) &(sorted_distances->front()),
			sizeof(double) * sorted_distances->size());
    nonsym_col += sorted_indices->size();      
    
  }


  // Final index -- just the total number of non-zero entries...
  nonsym_col_indices.write((char*) &nonsym_col, sizeof(int));

  // Close files
  distances.close();
  indices.close();
  nonsym_distances.close();
  nonsym_row_indices.close();
  nonsym_col_indices.close();

  o_file.write((char*) &current_index, (sizeof(int) / sizeof(char)));
  o_file.close();
  string mycall = "cat "
    + nsc_filename.str() + " "
    + nsr_filename.str() + " "
    + nsd_filename.str() + " >> " + o_filename;
  system(mycall.c_str());

  // Delete dynamically allocated data
  delete [] current_indices;
  delete [] current_distances;
  delete sorted_indices;
  delete sorted_distances;

  // Close the database
  if (cursor != NULL)
    cursor->close();
  db.close(0);

  // Delete the database file
  if (remove(db_file_name) != 0) {
    cout << "*** WARNING ***" << endl;
    cout << "Could not remove database file: " << db_file_name << endl;
  }
  if (remove(nsd_filename.str().c_str()) != 0) {
    cout << "*** WARNING ***" << endl;
    cout << "Could not remove database file: " << nsd_filename.str() << endl;
  }
  if (remove(nsr_filename.str().c_str()) != 0) {
    cout << "*** WARNING ***" << endl;
    cout << "Could not remove database file: " << nsr_filename.str() << endl;
  }
  if (remove(nsc_filename.str().c_str()) != 0) {
    cout << "*** WARNING ***" << endl;
    cout << "Could not remove database file: " << nsc_filename.str() << endl;
  }
  
  return 0;
}
