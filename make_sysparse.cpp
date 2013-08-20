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

// Berkeley DB
#include <db_cxx.h>

// Local
#include "config.h"

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
  
  if (argc != 2 && argc != 3) {
    cerr << endl;
    cerr << "   MDSCTK " << MDSCTK_VERSION_MAJOR << "." << MDSCTK_VERSION_MINOR << endl;
    cerr << "   Copyright (C) 2013 Joshua L. Phillips" << endl;
    cerr << "   MDSCTK comes with ABSOLUTELY NO WARRANTY; see LICENSE for details." << endl;
    cerr << "   This is free software, and you are welcome to redistribute it" << endl;
    cerr << "   under certain conditions; see README.md for details." << endl;
    cerr << endl;
    cerr << "Usage: " << argv[0] << " [k] <output k>" << endl;
    cerr << "   Converts the results from knn_rms into CSC format." << endl;
    cerr << endl;
    cerr << "   Normally, the number of nearest neighbors in the input" << endl;
    cerr << "   distances is used for constructing the CSC matrix." << endl;
    cerr << "   However, you can set <output k> <= [k] in order to" << endl;
    cerr << "   subselect the number of neighbors to consider in the" << endl;
    cerr << "   CSC representation. This makes it easy to store a" << endl;
    cerr << "   large number of neighbors using knn_* but then use" << endl;
    cerr << "   a subset for, say, computing approximate geodesic" << endl;
    cerr << "   distances." << endl;
    cerr << endl;
    return -1;
  }

  int k = 0;
  int maxk = 0;
  int current_index = 0;
  int *current_indices;
  double *current_distances;
  vector<int> *sorted_indices;
  vector<double> *sorted_distances;

  // Initialize data from command-line arguments
  maxk = atoi(argv[1]);
  if (argc == 4) {
    k = atoi(argv[2]);
  }
  else {
    k = maxk;
  }

  if (maxk < k) {
    cerr << "ERROR: Output k (" << k << ") is not less than the input k (" << maxk << ")." << endl;
    cerr << endl;
    return -1;
  }

  // char *distances_file_name = argv[3];
  // char *indices_file_name = argv[4];
  const char *distances_file_name = "distances.dat";
  const char *indices_file_name = "indices.dat";
  const char *sym_distances_file_name = "sym_distances.dat";
  const char *sym_row_indices_file_name = "sym_row_indices.dat";
  const char *sym_col_indices_file_name = "sym_col_indices.dat";

  int sym_col = 0;
  current_indices = new int[maxk];
  current_distances = new double[maxk];
  sorted_indices = new vector<int>;
  sorted_distances = new vector<double>;

  ifstream distances;
  ifstream indices;
  ofstream sym_distances;
  ofstream sym_row_indices;
  ofstream sym_col_indices;

  // Open files for reading/writing
  distances.open(distances_file_name,ios::in | ios::binary);
  if (!distances.good()) {
    cerr << "***ERROR***" << endl;
    cerr << "Could not open file: " << distances_file_name << endl;
    cerr << endl;
    return -1;
  }
  indices.open(indices_file_name,ios::in | ios::binary);
  if (!indices.good()) {
    cerr << "***ERROR***" << endl;
    cerr << "Could not open file: " << indices_file_name << endl;
    cerr << endl;
    distances.close();
    return -1;
  }
  sym_distances.open(sym_distances_file_name,ios::out | ios::binary);
  if (!sym_distances.good()) {
    cerr << "***ERROR***" << endl;
    cerr << "Could not open file: " << sym_distances_file_name << endl;
    cerr << endl;
    distances.close();
    indices.close();
    return -1;
  }
  sym_row_indices.open(sym_row_indices_file_name,ios::out | ios::binary);
  if (!sym_row_indices.good()) {
    cerr << "***ERROR***" << endl;
    cerr << "Could not open file: " << sym_row_indices_file_name << endl;
    cerr << endl;
    distances.close();
    indices.close();
    sym_distances.close();
    return -1;
  }
  sym_col_indices.open(sym_col_indices_file_name,ios::out | ios::binary);
  if (!sym_col_indices.good()) {
    cerr << "***ERROR***" << endl;
    cerr << "Could not open file: " << sym_col_indices_file_name << endl;
    cerr << endl;
    distances.close();
    indices.close();
    sym_distances.close();
    sym_row_indices.close();
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
  
  sprintf(db_file_name,"sym_%d.db",getpid());
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
    // cerr << "Starting work on " << current_index << ":";
    // for (int x = 0; x < k; x++)
    //   cerr << " " << current_indices[x];
    // cerr << endl;

    for (int x = 0; x < k; x++) {
            
      // Do Work
      myedge.from = current_index;	
      if (current_index < current_indices[x]) {
	myedge.to = current_indices[x];
	mydistance = current_distances[x];
	if (db.put(NULL, &key, &data, 0) != 0) {
	  cerr << "Could not insert edge: "
	       << myedge.from << " "
	       << myedge.to << " -> "
	       << mydistance << endl;
	}
      }	  

      myedge.to = current_index;
      if (current_indices[x] < current_index) {
	myedge.from = current_indices[x];
	mydistance = current_distances[x];
	if (db.put(NULL, &key, &data, 0) != 0) {
	  cerr << "Ccould not insert edge: "
	       << myedge.from << " "
	       << myedge.to << " -> "
	       << mydistance << endl;
	}
      }
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
    
    sym_col_indices.write((char*) &sym_col, (sizeof(int) / sizeof(char)));
    sym_row_indices.write((char*) &(sorted_indices->front()), 
			  (sizeof(int) / sizeof(char)) * sorted_indices->size());
    sym_distances.write((char*) &(sorted_distances->front()),
			(sizeof(double) / sizeof(char)) * sorted_distances->size());
    sym_col += sorted_indices->size();      
    
  }

  // Final index -- just the total number of non-zero entries...
  sym_col_indices.write((char*) &sym_col, (sizeof(int) / sizeof(char)));

  // Close files
  distances.close();
  indices.close();
  sym_distances.close();
  sym_row_indices.close();
  sym_col_indices.close();

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
    cerr << "*** WARNING ***" << endl;
    cerr << "Could not remove database file: " << db_file_name << endl;
  }
  
  return 0;
}
