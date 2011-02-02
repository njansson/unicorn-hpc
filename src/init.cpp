// Copyright (C) 2010 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2010-03-23
// Last changed: 2010-03-23


#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>
#include <dolfin.h>
#include "unicorn/unicorn_config.h"
#include "unicorn/UniParameters.h"
#include "unicorn/init.h"

using namespace dolfin;
using namespace dolfin::unicorn;

//-----------------------------------------------------------------------------
void unicorn::unicorn_init(int& argc, char* argv[], Mesh& mesh,
			   Checkpoint& chkp, long& w_limit, int& iter)
{
  dolfin_set("output destination","silent");
  if(dolfin::MPI::processNumber() == 0)
    dolfin_set("output destination","terminal");
  dolfin_init(argc, argv);
  message("Initializing Unicorn version %s.", PACKAGE_VERSION);

  if( argc < 5 ) 
    {
      message("Usage: -p <parameters> [-m <mesh> -c <checkpoint>] [-i iteration] [-l <wall clock limit>] [-o <petsc arguments>] ");      
      dolfin_finalize();
      exit(1);    
    }
  message("Running on %d %s", dolfin::MPI::numProcesses(), 
	  (dolfin::MPI::numProcesses() > 1 ? "nodes" : "node"));
  int c;
  while( -1 != (c = getopt(argc,argv,"p:m:c:i:l:o:")))
  {
    switch(c)
    {
    case 'p': 
      UniParameters::parse_parameters(optarg);
      break;
    case 'm':
      mesh = Mesh(optarg);
      break;
    case 'c':
      chkp.restart(optarg);  
      chkp.load(mesh);
      break;
    case 'i':
      iter = atoi(optarg);
      break;
    case 'o':
      break;
    case 'l':
      w_limit = atol(optarg);
      break;
    default:
      if(dolfin::MPI::processNumber() == 0)
        dolfin_set("output destination","terminal");      
      message("Unknown or missing argument");
      dolfin_finalize();
      exit(1);
      break;
    }
  }
  
  if(dolfin::MPI::processNumber() == 0)
    dolfin_set("output destination","terminal");  
  message("Global number of vertices: %d", 
	  (dolfin::MPI::numProcesses() > 1 ? mesh.distdata().global_numVertices() : mesh.numVertices()));
  message("Global number of cells: %d", 
	  (dolfin::MPI::numProcesses() > 1 ? mesh.distdata().global_numCells() : mesh.numCells()));
  dolfin_set("output destination","silent"); 
  
}
//-----------------------------------------------------------------------------
