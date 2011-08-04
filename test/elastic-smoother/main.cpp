
// Copyright (C) 2007 Johan Jansson.
// Licensed under the GNU GPL Version 2.

#include <dolfin.h>
#include <dolfin/mesh/RivaraRefinement.h>
#include <mpi.h>
#include "unicorn/ElasticSmoother.h"

using namespace dolfin;
using namespace dolfin::unicorn;

real bmarg = 1.0e-3 + DOLFIN_EPS;

namespace Geo
{
  
  // Geometry details /////////////////////////////////////////////////
  real cont_L = 3.0;
  real cont_H = 2.0;
  real cont_W = 2.0;
  
  real pin_l = 0.3;
  real pin_h = 0.005;
  real pin_w = 0.1;
  
  real circ_cx = 1.0;
  real circ_cy = 1.0;
  real circ_cz = 1.0;
  
  real pin_r = 0.05;
  
  //real pin_start_x = circ_cx;  
  real pin_start_x = circ_cx + pin_r;  
  real pin_end_x   = pin_start_x + pin_l;
  
  //real pin_start_y = circ_cy - (0.5 * pin_h);
  real pin_start_y = circ_cy + 0.045;
  real pin_end_y   = pin_start_y + pin_h;

  real pin_start_z = 0.95;
  real pin_end_z   = pin_start_z + pin_w;
  
  real xmin = 0.0; real xmax = cont_L;
  real ymin = 0.0; real ymax = cont_H;
  real zmin = 0.0; real zmax = cont_W;
  
  real pi = M_PI;
  
  real margin = DOLFIN_EPS + 1e-9;
  
  
  bool isStructure(Point r)
  {
    // Define the area that the pin occupies
    return (((pin_start_x - bmarg) < r[0] && r[0] < (pin_end_x + bmarg)) &&
	    ((pin_start_y - bmarg) < r[1] && r[1] < (pin_end_y + bmarg)) &&
	    ((pin_start_z - bmarg) < r[2] && r[2] < (pin_end_z + bmarg)));
  }
}


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

  Mesh mesh("mesh.xml");
  mesh.refine();
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

  for (CellIterator c(mesh); !c.end(); ++c)
  {
    Cell& cell = *c;
    Point mp = cell.midpoint();

    masked.set(cell, Geo::isStructure(mp));
  }

  MeshGeometry& geometry = mesh.geometry();

  for (VertexIterator v(mesh); !v.end(); ++v)
  {
    Vertex& vertex = *v;
    Point p = vertex.point();

    if(Geo::isStructure(p))
      geometry.x(vertex.index(), 1) += -0.003;
  }

  if(dolfin::MPI::processNumber() == 0)
    dolfin_set("output destination", "terminal");
  else
    dolfin_set("output destination", "silent");

  MeshQuality qual(mesh);
  qual.meshQuality();
  qual.disp();

  dolfin_set("Mesh read in serial", true);
  Mesh* structure_mesh = new Mesh("structure.xml");
  dolfin_set("Mesh read in serial", false);

  IntersectionDetector *idetector;
  idetector = new IntersectionDetector(*structure_mesh);

  MeshFunction<bool> solid_vertices(mesh, 0);

  for (VertexIterator v(mesh); !v.end(); ++v)
  {
    Vertex& vertex = *v;
    Point p = vertex.point();
    
    Array<unsigned int> overlap_cells;
    overlap_cells.clear();
    idetector->overlap(p, overlap_cells);
    
    bool bfnd = false;
    
    for(int i=0; i < overlap_cells.size();i++)
    {
      Cell testcell(*structure_mesh,overlap_cells[i]);
      if (structure_mesh->type().intersects(testcell,p))
      {
	std::cout << "solid vertex" << std::endl;
	bfnd = true;
	break;
      }			
    }
      
    solid_vertices.set(vertex, bfnd);
  }

  ElasticSmoother smoother(mesh);
  dolfin_set("Smoother max time steps", 20);
  smoother.smooth(smoothed, solid_vertices, h0);

  qual.meshQuality();
  qual.disp();

  return 0;
}
