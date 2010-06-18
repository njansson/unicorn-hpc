// Copyright (C) 2010 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2010-06-16
// Last changed: 2010-06-18

#ifndef __DATASETSAMPLER_H
#define __DATASETSAMPLER_H

#include <dolfin.h>

namespace dolfin
{
  namespace unicorn 
  {
    class DataSetSampler
    {
    public:
      
      // Sample dataset on a cartesian grid
      static void cartesian_sampling(Function& f, real *value, uint N,
				     real xmin, real xmax,
				     real ymin, real ymax,
				     real zmin, real zmax);
    };
  }
}
#endif
