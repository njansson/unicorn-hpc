// Copyright (C) 2005 Johan Hoffman.
// Licensed under the GNU GPL Version 2.
//
// Modified by Anders Logg 2005.
// Modified by Niclas Jansson 2008-2010.
//
// First added:  2002-11-29
// Last changed: 2010-06-13
//
// A cG(1)cG(1) FEM solver for the incompressible Navier-Stokes equations 
//
//     du/dt + u * grad u - nu * div grad u + grad p = f 
//     div u = 0 
//

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>

#include <dolfin.h>
#include <unicorn/init.h>
#include <unicorn/util.h>
#include <unicorn/NodeNormal.h>
#include <unicorn/UniParameters.h>
#include <unicorn/NSESolver.h>
#include <unicorn/NewSlipBC.h>
#include <unicorn/NewDirichletBC.h>

using namespace dolfin;
using namespace dolfin::unicorn;

real bmarg = 1.0e-3 + DOLFIN_EPS;

namespace Geo
{
  
  // Geometry details
  real cont_L = 2.5;
  real cont_H = 0.41;
  
  real pin_l = 0.35;
  real pin_h = 0.02;
  
  real circ_cx = 0.2;
  real circ_cy = 0.2;
  
  real pin_r = 0.05;
  
  real pin_start_x = circ_cx;  
  real pin_end_x   = pin_start_x + pin_r + pin_l;
  
  real pin_start_y = circ_cy - (0.5 * pin_h);
  real pin_end_y   = pin_start_y + pin_h;
  
  real xmin = 0.0; real xmax = cont_L;
  real ymin = 0.0; real ymax = cont_H;
  
  real pi = M_PI;
  
  real margin = DOLFIN_EPS + 1e-9;
  
  
  bool isStructure(Point r)
  {
    // // Define the area that the pin occupies
    // // return (((pin_start_x - margin) < r[0] && r[0] < (pin_end_x + margin)) &&
    // // 	    ((pin_start_y - margin) < r[1] && r[1] < (pin_end_y + margin)));
    // real val =
    //   (((pin_start_x - margin) < r[0] && r[0] < (pin_end_x + margin)) &&
    //    ((pin_start_y - margin) < r[1] && r[1] < (pin_end_y + margin)));
    // std::cout << "isStructure: " << val << std::endl;
    // return val;
    return true;
  }
}

real xmin = 0.0;
real xmax = 2.2;
//real xmax = 1.0;
real ymin = 0.0;
real ymax = 0.41;
//real ymax = 1.0;
real zmin = 0.0;
real zmax = 0.4;
real xcenter = 0.2;
real ycenter = 0.2;
real radius = 0.05;

real density(Point p)
{
  real val = 0.0;

  val = 1.1 - 0.1 * p[1];

  /*
  if(p[1] > 0.5)
    {
      if(2.0 * p[0] > p[1])
	{
	  val = 1.1;
	}
      else
	{
	  val = 1.0;
	}
    }
  else
    {
      val = 1.1;
    }
  */
  return val;
}

class AllBoundary : public SubDomain
{
public:
  bool inside(const real* p, bool on_boundary) const
  {
    return on_boundary;
  }
};

// Sub domain for inflow
class InflowTopBottomBoundary3D : public SubDomain
{
public:
  bool inside(const real* p, bool on_boundary) const
  {
    return on_boundary && (p[0] < xmax - bmarg ||
                           p[1] < ymin + bmarg || p[1] > ymax - bmarg);
  }
};

// Sub domain for outflow
class OutflowBoundary : public SubDomain
{
public:
  bool inside(const real* p, bool on_boundary) const
  {
    InflowTopBottomBoundary3D ib;

    return on_boundary && !ib.inside(p, on_boundary);
  }
};

// Sub domain for slip
class SlipBoundary : public SubDomain
{
public:
  bool inside(const real* p, bool on_boundary) const
  {

    return on_boundary && false && p[0] < (xmax - bmarg) && p[0] > (xmin + bmarg);
  }
};

// Sub domain for Channel side walls
class SideWallBoundary : public SubDomain 
{ 
public: 
  bool inside(const real* p, bool on_boundary) const 
  { 
    return on_boundary && 
      (p[0] > xmin + bmarg && p[0] < xmax - bmarg) && 
      ((p[1] < ymin + bmarg || p[1] > ymax - bmarg) || 
       (p[2] < zmin + bmarg || p[2] > zmax - bmarg)); 
  } 
}; 

// Force term
class ForceFunction : public Function
{
public:
  ForceFunction(Mesh& mesh) : Function(mesh) {}
  void eval(real* values, const real* x) const
  {
    int d = cell().dim();

    for(int i = 0; i < d; i++)
    {
      values[i] = 0.0;
    }
  }
};

// Boundary condition for momentum equation
class BC_Momentum_3D : public Function
{
public:
  BC_Momentum_3D(Mesh& mesh, TimeDependent& td) : Function(mesh), td(td) {}
  void eval(real* values, const real* x) const
  {
    real t = td.time();
    real ramp = 0.0;
    real tramp = 1.0;

    if(t < tramp)
    {
      ramp = (t - 0.0) / (tramp - 0.0);
    }
    else if(t >= tramp)
    {
      ramp = 1.0;
    }

    for(int i = 0; i < 2; i++)
    {
      if (i==0)
      {
	if ( x[0] < (0.0 + bmarg) &&
	     x[1] < (ymax - bmarg) &&
	     x[1] > (ymin + bmarg))
	  {
	  real Um = 2.0;
	  
	  
	  values[i] =
	    ramp*1.5*Um*( (x[1]*(Geo::cont_H - x[1])) /
			  ((0.5*Geo::cont_H)*(0.5*Geo::cont_H)) );
	  }
	else
	{
	  values[i] = 0.0;
	}
      } 
      else if (i==1)
      {
	values[i] = 0.0;
      } 
    }
  }

  TimeDependent& td;
};

// Boundary condition for dual momentum equation
class Dual_BC_Momentum_3D : public Function
{
public:
  Dual_BC_Momentum_3D(Mesh& mesh) : Function(mesh) {}
  void eval(real* values, const real* x) const
  {
    if ( (x[0] > xmin + bmarg) && (x[0] < xmax - bmarg) &&
	 (x[1] > ymin + bmarg) && (x[1] < ymax - bmarg) &&
	 (x[2] > zmin + bmarg) && (x[2] < zmax - bmarg) ) {
      values[0] = -1.0;
      values[1] = 0.0;
      values[2] = 0.0;
    }
    else {
      values[0] = 0.0;
      values[1] = 0.0;
      values[2] = 0.0;
    }      
  }
};


// Boundary condition for continuity equation 
class BC_Continuity : public Function
{
public:
  BC_Continuity(Mesh& mesh) : Function(mesh) {}
  
  void eval(real* values, const real* x) const
  {
    values[0] = 0.0;
  }
};


class Phi : public Function
{
public:
  Phi(Mesh& mesh) : Function(mesh) {}
  void eval(real* values, const real* x) const
  {
    int d = cell().dim();

    for(int i = 0; i < d; i++)
    {
      values[i] = 0.0;
    }

    if ( (x[0] > xmin + bmarg) && (x[0] < xmax - bmarg) &&
	 (x[1] > ymin + bmarg) && (x[1] < ymax - bmarg) &&
	 (x[2] > zmin + bmarg) && (x[2] < zmax - bmarg) )
      values[0] = 1.0;
  }
};

// Friction coefficient
class Beta : public Function
{
public:
  Beta(Mesh& mesh) : Function(mesh) {}
  void eval(real* values, const real* x) const
  {
    if(x[0] > xmin + bmarg && x[0] < xmax - bmarg &&
       x[1] > ymin + bmarg && x[1] < ymax - bmarg &&
       x[2] > zmin + bmarg && x[2] < zmax - bmarg)
    {
      values[0] = dolfin_get("beta");
    }
    else
    {
      values[0] = 0.0;
    }
  }
};

class BCDensity : public Function
{
public:
  BCDensity(Mesh& mesh) : Function(mesh) {}
  void eval(real* values, const real* x) const
  {
    Point p(x[0], x[1], x[2]);

    values[0] = density(p);
  }
};


void solve(Mesh& mesh, Checkpoint& chkp, long& w_limit, timeval& s_time, Mesh* structure_mesh)
{
  real T = dolfin_get("T");
  real dual_T = dolfin_get("dual_T");
  real nu = dolfin_get("nu");
  real ubar = dolfin_get("Ubar");
  real mu = dolfin_get("mu");
  real nu_s = dolfin_get("nu_s");
  
  TimeDependent td;

  InflowTopBottomBoundary3D iboundary;
  BC_Momentum_3D bcf_mom(mesh, td);
  BCDensity density_function(mesh);
  Dual_BC_Momentum_3D dual_bcf_mom(mesh);
  BC_Continuity bcf_con(mesh);
  ForceFunction f(mesh);
  Function zero(mesh, 2, 0.0);
  AllBoundary aboundary;
  OutflowBoundary oboundary;
  SlipBoundary slipboundary;
  SideWallBoundary swboundary;

  Phi phi(mesh);
  Beta beta(mesh);

  Array<Function*> aero_f;
  aero_f.push_back(&phi);

  //NodeNormal node_normal(mesh);
  DirichletBC p_bcin(bcf_mom, mesh, iboundary);
  DirichletBC p_bcout(bcf_con, mesh, oboundary);
  DirichletBC p_bcdensity(density_function, mesh, iboundary);
  //NewSlipBC slip_bc(mesh, slipboundary, node_normal);
  
  //dolfin_set("PDE slip alpha max", DOLFIN_PI / 1.9);
  Array <BoundaryCondition*> p_bc_momentum;
  p_bc_momentum.push_back(&p_bcin);
  //p_bc_momentum.push_back(&slip_bc);

  Array <BoundaryCondition*> p_bc_pressure;
  p_bc_pressure.push_back(&p_bcout);
  Array <BoundaryCondition*> p_bc_density;
  p_bc_density.push_back(&p_bcdensity);
  //p_bc_density.push_back(&p_bcout);

  MeshFunction<bool> solid_cells(mesh, mesh.topology().dim());

  bool existRegionFile = structure_mesh;
   
  IntersectionDetector *idetector;
   
  if (existRegionFile) 
  {
    std::cout << "Region file exists" << std::endl;
    idetector = new IntersectionDetector(*structure_mesh);
  }
  else
  {
    std::cout << "Region file doesn't exist" << std::endl;
  }

  for (CellIterator c(mesh); !c.end(); ++c)
  {
    Cell& cell = *c;
    Point mp = cell.midpoint();

   if (existRegionFile)
   {
     Array<unsigned int> overlap_cells;
     overlap_cells.clear();
     idetector->overlap(mp, overlap_cells);
     
     bool bfnd = false;
     for(int i=0; i < overlap_cells.size(); i++)
     {	
       Cell testcell(*structure_mesh, overlap_cells[i]);
       if (structure_mesh->type().intersects(testcell,mp))
       {
	 std::cout << "solid cell" << std::endl;
	 bfnd = true;
	 break;
       }			
     }
     
     solid_cells.set(cell, bfnd);
     
   }
   else
     solid_cells.set(cell, Geo::isStructure(mp));
   //solid_cells.set(cell, false);
   if(solid_cells.get(cell))
   {
     std::cout << "solid0" << std::endl;
   }
  }
  
  MeshFunction<bool> solid_vertices(mesh, 0);
  
  for (VertexIterator v(mesh); !v.end(); ++v)
  {
    Vertex& vertex = *v;
    Point p = vertex.point();
    
    if (existRegionFile)
    {
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
    else
      solid_vertices.set(vertex, Geo::isStructure(p));
  }

  // MeshFunction<int> old2new_vertex;
  // MeshFunction<int> old2new_cells;
  // Mesh sub;
  // ElasticSmoother::submesh(mesh, sub, solid_cells, old2new_vertex, old2new_cells);

  // File submeshfile("sub.xml");
  // File submeshfile2("sub.pvd");
  // submeshfile << sub;
  // submeshfile2 << sub;

  real rho_s = 1.0;

  Function U, U0;

//   dolfin_set("Krylov relative tolerance", 1.0e-12);
//   dolfin_set("Krylov absolute tolerance", 1.0e-20);

  NSESolver psolver(mesh, U, U0, f, f, phi, beta, p_bc_momentum, p_bc_pressure,
		    p_bc_density, &density, solid_cells, solid_vertices, T, nu, mu, nu_s, rho_s,
                    ubar, td, "primal"); 
  dolfin_set("Adaptive refinement percentage", 5.0);
  dolfin_set("ODE discrete tolerance", 1.0e-2);
  dolfin_set("ODE maximum iterations", 30);
  dolfin_set("PDE number of samples", 500);  

  psolver.solve(U, U0);
  
}


int main(int argc, char* argv[])
{
  timeval s_time;
  gettimeofday(&s_time, NULL);
  Mesh mesh;  
  long w_limit = 0;
  Checkpoint chkp;
  int iter = 0;
  Mesh* structure_mesh;

  unicorn_init(argc, argv, mesh, chkp, w_limit, iter, structure_mesh);

  //mesh.refine();
  //mesh.refine();

  unicorn_solve(mesh, chkp, w_limit, s_time, iter, 0, 0, &solve, structure_mesh);

  dolfin_finalize();
  return 0;
}
