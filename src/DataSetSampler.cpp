// Copyright (C) 2010 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2010-06-16
// Last changed: 2010-06-18

#include "unicorn/DataSetSampler.h"

using namespace dolfin;
using namespace unicorn;

//-----------------------------------------------------------------------------
void DataSetSampler::cartesian_sampling(Function& f, real *value, uint N,
					real xmin, real xmax,
					real ymin, real ymax,
					real zmin, real zmax)
{
  
  real x[3];
  real origin[3];
  real space = (xmax - xmin) / N;
  origin[0] = xmin;
  origin[1] = ymin;
  origin[2] = zmin;

  real *vp = &value[0];
  for (uint k = 0; k < N; k++)  
  {       	
    x[2] = k*space + zmin;
    for (uint j = 0; j < N; j++)  
    {
      x[1] = j*space + ymin;
      for (uint i = 0; i < N; i++)  
      {
	x[0] = i*space + xmin;
	f.eval((vp++), &x[0]);
      }
    }
  }
  
#ifdef UNICORN_DEBUG
  FILE *out;
  out = fopen("test.vtk","w");
  fprintf(out,"# vtk DataFile Version 3.0\n");
  fprintf(out,"Sampling\n");
  fprintf(out,"ASCII\n");
  fprintf(out,"DATASET STRUCTURED_POINTS\n");
  fprintf(out,"DIMENSIONS %d %d %d\n",N,N,N);
  fprintf(out,"ORIGIN %e %e %e\n",origin[0],origin[1],origin[2]);
  fprintf(out,"SPACING %e %e %e\n",space,space,space);
  fprintf(out,"POINT_DATA %d\n",N*N*N);
  fprintf(out,"SCALARS Sampling double 1\n");
  fprintf(out,"LOOKUP_TABLE default\n");
  
  for(int i=0; i<N*N*N; i++) {
    fprintf(out,"%e\n",value[i]);
  }
  fclose(out);
#endif
}
//-----------------------------------------------------------------------------
