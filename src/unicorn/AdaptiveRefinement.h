// Copyright (C) 2010 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2010-09-13
// Last changed: 2010-10-11

#ifndef __ADAPTIVEREFINEMENT_H
#define __ADAPTIVEREFINEMENT_H

#include <vector>
#include <dolfin.h>

namespace dolfin
{
  namespace unicorn 
  {
    class AdaptiveRefinement
    {
    public:
      

      typedef std::pair<Form *, uint> form_tuple;
      typedef std::pair< Function *, form_tuple > project_func;
      
      static void refine(Mesh& mesh, MeshFunction<bool>& cell_marker);

      static void refine_and_project(Mesh& mesh, std::vector<project_func> pf,
				     MeshFunction<bool>& cell_marker);

    private:

      static void redistribute_func(Mesh& mesh, Function *f, 
				    real **vp, uint **rp, uint& m,
				    Form& form, uint offset,
				    MeshFunction<uint>& distribution);

      static void decompose_func(Mesh& mesh, Function *f, uint offset, Form& form,
				 Function& f_x, Function& f_y, Function& f_z);


      // Comparison operator for index/value pairs    
      struct less_pair : public std::binary_function<std::pair<uint, real>,
						     std::pair<uint, real>, bool>
      {
        bool operator()(std::pair<uint, real> x, std::pair<uint, real> y)
	{
	  return x.second < y.second;
	}
      };
    

      
    };
  }
}
#endif
