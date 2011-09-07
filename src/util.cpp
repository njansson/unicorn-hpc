// Copyright (C) 2010 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2010-03-29
// Last changed: 2010-03-29

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <dolfin.h>
#include "unicorn/UniParameters.h"
#include "unicorn/util.h"

using namespace dolfin;
using namespace dolfin::unicorn;

//-----------------------------------------------------------------------------
void unicorn::unicorn_solve(Mesh& mesh, Checkpoint& chkp, 
			    long& w_limit, timeval& s_time, int& iter,
			    void (*pre)(Mesh&), void (*post)(Mesh&),
			    void (*solver)(Mesh&, Checkpoint&, long&, timeval&, Mesh*),
			    Mesh* structure_mesh)
{
  char itername[80];
  for(int i = iter; i < (int) dolfin_get("adapt_iter") ; i++) 
    {
      if(dolfin::MPI::processNumber() == 0)
	dolfin_set("output destination", "terminal");
      message("Running iteration %d of %d", i, (int) dolfin_get("adapt_iter"));
      dolfin_set("output destination", "silent");
      snprintf(itername, sizeof(itername), "iter_%d", i);
      
      if (dolfin::MPI::processNumber() == 0) 
	if(mkdir(itername, S_IRWXU) < 0)
	  perror("mkdir failed");
      
      MPI_Barrier(dolfin::MPI::DOLFIN_COMM);
      
      if(chdir(itername) < 0)
	perror("chdir failed");
      
      if (pre)
      	(*pre)(mesh);
	  
      (*solver)(mesh, chkp, w_limit, s_time, structure_mesh);

      if (post)
      	(*post)(mesh);
      
      if(chdir("../") < 0)
	perror("chdir failed");
      
    }  
}
//-----------------------------------------------------------------------------
