
// Copyright (C) 2007 Johan Jansson.
// Licensed under the GNU GPL Version 2.

#include <dolfin.h>
#include <dolfin/mesh/RivaraRefinement.h>
#include <mpi.h>
#include "unicorn/ElasticSmoother.h"

using namespace dolfin;
using namespace dolfin::unicorn;

void transform(Mesh& mesh)
{
  for (VertexIterator n(mesh); !n.end(); ++n)
  {
    Vertex& vertex = *n;
    Point p = vertex.point();
      
    MeshGeometry& geometry = mesh.geometry();
      
    geometry.x(vertex.index(), 0) *= 0.01;
    geometry.x(vertex.index(), 1) *= 0.01;
    geometry.x(vertex.index(), 2) *= 0.01;
  }
}


int main(int argc, char* argv[])
{
  dolfin_init(argc, argv);

  //UnitSquare mesh(20, 2);
  //UnitCube mesh(20, 2, 2);
  
  //Mesh mesh("testcell3D.xml");
  //Mesh mesh("cubemesh.xml");
  //Mesh mesh("test.xml");
  //UnitCube mesh(3, 3, 3);
  //Mesh mesh("heart_coarse.xml");
  //Mesh mesh("heart_fine.xml");
  //Mesh mesh("heart_new03.xml");
  //Mesh mesh("tire.xml");
  Mesh mesh("wedge.xml");

  //transform(mesh);

  mesh.renumber();

  for(int i = 0; i < 0; i++)
  {
    
    MeshFunction<bool> cell_refinement_marker(mesh);
    cell_refinement_marker.init(mesh.topology().dim());
    cell_refinement_marker = false;
    
    for (CellIterator c(mesh); !c.end(); ++c)
    {
      if(c->midpoint()[0] < 0.25)
	//if(c->midpoint()[0] > 0.75)
	cell_refinement_marker.set(c->index(), true);
    }
    
    RivaraRefinement::refine(mesh, cell_refinement_marker);
  }
  //Mesh mesh("testcell3D.xml");

  MeshFunction<bool> smoothed(mesh, mesh.topology().dim());
  MeshFunction<bool> masked(mesh, mesh.topology().dim());
  MeshFunction<real> h0(mesh, mesh.topology().dim());

  smoothed = true;
  masked = false;
  
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    Cell& cell = *c;

    //h0.set(cell, cell.diameter());
    h0.set(cell, MeshQuality::myDiameter(cell));
  }

  ElasticSmoother smoother(mesh);
  smoother.smooth(smoothed, masked, h0);

  return 0;
}
