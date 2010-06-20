// Copyright (C) 2007 Murtazo Nazarov
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2009.
// 
// First added:  2007-05-01
// Last changed: 2009-03-17
                                                                                                  
#ifndef __NODENORMAL_H
#define __NODENORMAL_H
#include <dolfin/mesh/Mesh.h>

#include <dolfin/common/constants.h>
#include <dolfin/common/Array.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <map>

namespace dolfin
{
  namespace unicorn
  {
    class NodeNormal
    {
    public:
      
      // Copy constructor
      NodeNormal(NodeNormal& node_normal);
      
      // Create normal, tangents for the boundary of mesh
      NodeNormal(Mesh& mesh);
      
      ~NodeNormal();
      
      // Assignment
      NodeNormal& operator=(NodeNormal& node_normal);
      // Cleanup
      void clear();
      
      // Compute normals to the boundary nodes
      void ComputeNormal(Mesh& mesh);
      
      Mesh& mesh;    
      
      // Define mesh functions for normal and tangents
      MeshFunction<real> *normal, *tau, *tau_1, *tau_2;
      
      // Define node type: 1 surface, 2 edge, 3 surface
      MeshFunction<int> node_type;
      
    private:
      
      void cache_shared_area(Mesh& mesh, BoundaryMesh& boundary, uint nsdim, 
			     MeshFunction<uint> *vertex_map, 
			     MeshFunction<uint> *cell_map);
      
      std::map<uint, uint> shared_noffset;
      std::map<uint, uint> num_cells;
      Array<real> shared_normal;    
      std::map<uint, Array<real> > normal_block;
      std::map<uint, Array<real> > shared_area_block;
      uint _offset;
      
    };    
  }
}
#endif
  
