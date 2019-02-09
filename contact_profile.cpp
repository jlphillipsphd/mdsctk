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

int main(int argc, char* argv[]) {

  const char* program_name = "contact_profile";
  bool optsOK = true;
  gmx::initForCommandLine(&argc,&argv);
  copyright(program_name);
  cout << "   Computes the standard atomic contacts for structures in" << endl;
  cout << "   the given xtc file. A topology PDB file and atom index file" << endl;
  cout << "   should be provided for determining the atoms to compare." << endl;
  cout << "   The resulting sparse contact distance profiles are" << endl;
  cout << "   in sparse vector format (index-file and data-file)." << endl;
  cout << endl;
  cout << "   Use -h or --help to see the complete list of options." << endl;
  cout << endl;

  // Option vars...
  int nthreads = 0;
  double sigma;
  double eps;
  string top_filename;
  string xtc_filename;
  string ndx_filename;
  const char* ndx_filename_ptr = NULL;
  string index_filename;
  string data_filename;

  // Declare the supported options.
  po::options_description cmdline_options;
  po::options_description program_options("Program options");
  program_options.add_options()
    ("help,h", "show this help message and exit")
    ("threads,t", po::value<int>(&nthreads)->default_value(omp_get_max_threads()>omp_get_num_procs()?omp_get_num_procs():omp_get_max_threads()), "Input:  Number of threads to start (int)")
    ("epsilon,e", po::value<double>(&eps)->default_value(9.0), "Input:  Contact cutoff (real)")
    //    ("sigma,q", po::value<double>(&sigma)->default_value(1), "Input:  Standard deviation of gaussian kernel (real)")
    ("topology-file,p", po::value<string>(&top_filename)->default_value("topology.pdb"), "Input:  Topology file [.pdb,.gro,.tpr] (string:filename)")
    ("xtc-file,x", po::value<string>(&xtc_filename)->default_value("traj.xtc"), "Input:  Trajectory file (string:filename)")
    ("ndx-file,n", po::value<string>(&ndx_filename), "Input: K-nn distances file (string:filename)")
    ("index-file,i", po::value<string>(&index_filename)->default_value("reference.svi"), "Output: Sparse vector indices file (string:filename)")    
    ("data-file,d", po::value<string>(&data_filename)->default_value("reference.svd"), "Output: Sparse vector data file (string:filename)")    
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
  if (vm.count("ndx-file")) {
    ndx_filename_ptr = ndx_filename.c_str();
  }

  if (!optsOK) {
    return -1;
  }

  cout << "Running with the following options:" << endl;
  cout << "threads =       " << nthreads << endl;
  cout << "topology-file = " << top_filename << endl;
  cout << "xtc-file =      " << xtc_filename << endl;
  cout << "ndx-file =      " << ndx_filename << endl;
  cout << "index-file =    " << index_filename << endl;
  cout << "data-file =     " << data_filename << endl;
  cout << endl;
  
  // Local vars
  long int step = 1;
  float time = 0.0;
  matrix box;
  float prec = 0.001;
  // char buf[256];
  t_topology top;
  int ePBC;
  int natoms = 0;
  int nframes= 0;
  int update_interval = 1;
  t_fileio *ref_file;
  rvec *mycoords = NULL;
  gmx_bool bOK = 1;
  double *contact = NULL;
  vector<coord_array> *ref_coords = NULL;
  ::real *weights = NULL;
  int        gnx1,gnx2;
  int        *index1,*index2;
  char       *grpname1,*grpname2;
  ofstream   index;
  ofstream   data;

  // Remove C stdout (silly GROMACS warnings going every which stream!)
  int myout = dup(1);
  dup2(2,1);

  // Setup threads
  omp_set_num_threads(nthreads);

  // Get number of atoms and check xtc
  cout << "Reading topology information from " << top_filename << " ... ";
  read_tps_conf(top_filename.c_str(), &top, &ePBC, &mycoords,
		NULL, box, TRUE);
  cout << "done." << endl;
  delete [] mycoords;

  ref_file = open_xtc(xtc_filename.c_str(),"r");
  read_first_xtc(ref_file,&natoms, &step, &time, box, &mycoords, &prec, &bOK);
  close_xtc(ref_file);
  if (natoms != top.atoms.nr) {
    cout << "*** ERROR ***" << endl;
    cout << "Number of atoms in topology file ("
	 << top.atoms.nr << ") "
	 << "does not match the number of atoms "
	 << "in the XTC file (" << xtc_filename << " : " << natoms << ")."
	 << endl;
    exit(4);
  }

  // Get atom selections
  cout << "Please select two (non-overlapping) groups for contact profiling..." << endl;
  get_index(&top.atoms,ndx_filename_ptr,1,&gnx1,&index1,&grpname1);
  cout << endl;
  get_index(&top.atoms,ndx_filename_ptr,1,&gnx2,&index2,&grpname2);
  cout << endl;

  cout << "Total grid size is " << gnx1 << " x " << gnx2 << " = " << (gnx1*gnx2) << endl;

  // Read coordinates and weight-center all structures
  cout << "Reading reference coordinates from file: " << xtc_filename << " ... ";
  ref_coords = new vector<coord_array>;
  ref_file = open_xtc(xtc_filename.c_str(),"r");
  mycoords = new rvec[natoms];
  while (read_next_xtc(ref_file, natoms, &step, &time, box, mycoords, &prec, &bOK)) {
    ref_coords->push_back(mycoords);
    mycoords = new rvec[natoms];
  }
  close_xtc(ref_file);
  delete [] mycoords;
  mycoords = NULL;
  nframes = ref_coords->size();
  cout << "done." << endl;

  // Allocate vectors for storing the distances for a structure
  contact = new double[gnx1*gnx2];
  weights = new ::real[gnx1*gnx2];
  for (int x = 0; x < natoms; x++) weights[x] = top.atoms.atom[x].m;

#pragma omp parallel for
  for (int i = 0; i < gnx1; i++)
    for (int j = 0; j < gnx2; j++) {
      weights[(i*gnx2)+j] = top.atoms.atom[index1[i]].m * top.atoms.atom[index2[j]].m;
    }

  // Restore C stdout.
  dup2(myout,1);

  index.open(index_filename.c_str());
  data.open(data_filename.c_str());

  // Timer for ETA
  time_t start = std::time(0);
  time_t last = start;

  // Compute fits
  for (int frame = 0; frame < nframes; frame++) {

    // Update user of progress
    if (std::time(0) - last > update_interval) {
      last = std::time(0);
      time_t eta = start + ((last-start) * nframes / frame);
      cout << "\rFrame: " << frame << ", will finish " 
	   << string(std::ctime(&eta)).substr(0,20);
      cout.flush();
    }

    // Do Work
#pragma omp parallel for
    for (int i = 0; i < gnx1*gnx2; i++)
      contact[i] = 0.0;

#pragma omp parallel for
    for (int i = 0; i < gnx1; i++) {
      int ii = index1[i];
      for (int j = 0; j < gnx2; j++) {
	int jj = index2[j];
	double d = 0.0;
	for (int k = 0; k < 3; k++)
	  d += (((*ref_coords)[frame][ii][k] - (*ref_coords)[frame][jj][k]) *
		((*ref_coords)[frame][ii][k] - (*ref_coords)[frame][jj][k]));
	d = sqrt(d) * 10.0;
	// d = exp(-(d*d) / (2.0 * weights[(i*gnx2)+j]));
	// if (d > eps)
	//   contact[(i*gnx2)+j] = d;
	if (d < eps)
	  contact[(i*gnx2)+j] = 1.0;
	else if (eps <= 0.0)
	  contact[(i*gnx2)+j] += d;
      } // j
    } // i


    double sum = 0.0;
#pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < gnx1*gnx2; i++)
      sum += contact[i];

    sum = 1.0; // No normalization...
    int total = 0;
#pragma omp parallel for reduction(+:total)
    for (int i = 0; i < gnx1*gnx2; i++)
      if (contact[i] > 0) {
	contact[i] /= sum;
	total++;
      }

    index.write((char*) &total, sizeof(int) / sizeof(char));
    for (int i = 0; i < gnx1*gnx2; i++)
      if (contact[i] > 0.0) {
	index.write((char*) &i, sizeof(int) / sizeof(char));
	data.write((char*) &contact[i], sizeof(double) / sizeof(char));
      }

    // cout << frame << " " << total << endl;

  } // frame

  cout << endl << endl;

  index.close();
  data.close();

  // Clean coordinates
  for (vector<coord_array>::iterator itr = ref_coords->begin();
       itr != ref_coords->end(); itr++) delete [] (*itr);
  delete ref_coords;

  delete [] contact;
  delete [] weights;
  
  return 0;
}
