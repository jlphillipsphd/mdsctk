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

int main(int argc, char *argv[]) {

  const char* program_name = "make_gesparse";
  bool optsOK = true;
  gmx::initForCommandLine(&argc,&argv);
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
  bool s;
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
    ("symmetric,s", po::bool_switch(&s)->default_value(false), "Input:  Enforce symmetry (bool)")
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
    cout << "Could not open file: " << nsd_filename.str() << endl;
    cout << endl;
    return -1;
  }
  nonsym_row_indices.open(nsr_filename.str().c_str(),ios::out | ios::binary);
  if (!nonsym_row_indices.good()) {
    cout << "***ERROR***" << endl;
    cout << "Could not open file: " << nsr_filename.str() << endl;
    cout << endl;
    return -1;
  }
  nonsym_col_indices.open(nsc_filename.str().c_str(),ios::out | ios::binary);
  if (!nonsym_col_indices.good()) {
    cout << "***ERROR***" << endl;
    cout << "Could not open file: " << nsc_filename.str() << endl;
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

  cout << "Creating sparse matrix database..." << endl;

  // Timer for ETA
  time_t start = std::time(0);
  time_t last = start;
  int nframes = 0;
  indices.seekg(0,ios::end);
  nframes = (indices.tellg() * sizeof(char) / sizeof(int)) / maxk;
  indices.seekg(0,ios::beg);

  // Read a set of distances and indices
  distances.read((char*) current_distances, sizeof(double) * maxk);
  indices.read((char*) current_indices, sizeof(int) * maxk);

  while (!distances.eof() || !indices.eof()) {
 
    // Update user of progress
    if (std::time(0) - last > update_interval) {
      last = std::time(0);
      time_t eta = start + ((last-start) * nframes / current_index);
      cout << "\rFrame: " << current_index << ", will finish " 
	   << string(std::ctime(&eta)).substr(0,20);
      cout.flush();
    }
    
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
      if (s) {
	myedge.to = current_index;
	myedge.from = current_indices[x];
	mydistance = current_distances[x];
	if (db.get(NULL, &key, &data, 0) == DB_NOTFOUND) {
	  myedge.to = current_index;
	  myedge.from = current_indices[x];
	  mydistance = current_distances[x];
	  if (db.put(NULL, &key, &data, 0) != 0) {
	    cout << "Ccould not insert edge: "
		 << myedge.from << " "
		 << myedge.to << " -> "
		 << mydistance << endl;
	  }
	}
      }
    }

    // Read a set of distances and indices
    distances.read((char*) current_distances, sizeof(double) * maxk);
    indices.read((char*) current_indices, sizeof(int) * maxk);
    current_index++;
  }

  cout << endl;
  cout << "Converting database to sparse matrix..." << endl;

  start = std::time(0);
  last = start;
  nframes = current_index;

  // Initialize cursor
  db.cursor(NULL, &cursor, 0); 
  cursor->get(&key, &data, DB_FIRST);
    
  for (int col = 0; col < current_index; col++) {

    // Update user of progress
    if (std::time(0) - last > update_interval) {
      last = std::time(0);
      time_t eta = start + ((last-start) * nframes / col);
      cout << "\rFrame: " << col << ", will finish " 
	   << string(std::ctime(&eta)).substr(0,20);
      cout.flush();
    }
    
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
  if (system(mycall.c_str()))
    cout << "Could not create general CSC matrix file: " << o_filename << endl;

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
