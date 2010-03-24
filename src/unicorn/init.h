// Copyright (C) 2010 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2010-03-23
// Last changed: 2010-03-23

#ifndef __UNICORN_INIT_H
#define __UNICORN_INIT_H

#define UNICORN_VERSION PACKAGE_VERSION

#include <dolfin.h>

namespace dolfin
{
  namespace unicorn
  {
    
    void unicorn_init(int& argc, char* argv[], Mesh& mesh, 
		      Checkpoint& chkp, long& w_limit, int& iter);
  }
}
#endif
