// Copyright (C) 2010 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2010-03-29
// Last changed: 2010-03-29

#ifndef __UNICORN_UTIL_H
#define __UNICORN_UTIL_H

#include <sys/time.h>
#include <dolfin.h>

namespace dolfin
{
  namespace unicorn
  {
    
    void unicorn_solve(Mesh& mesh, Checkpoint& chkp, 
		       long& w_limit, timeval& s_time, int& iter,
		       void (*pre)(Mesh&), void (*post)(Mesh&),
		       void (*solver)(Mesh&, Checkpoint&, long&, timeval&, Mesh*),
		       Mesh* structure_mesh = 0);
  }
}
#endif
