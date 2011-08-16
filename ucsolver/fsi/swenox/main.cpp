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
#include <unicorn/SlipBC.h>
#include <unicorn/NewDirichletBC.h>

using namespace dolfin;
using namespace dolfin::unicorn;

real bmarg = 1.0e-5 + DOLFIN_EPS;
real gmarg = 1.0e-5 + DOLFIN_EPS;

namespace Geo
{
  
  // Geometry details /////////////////////////////////////////////////
  real cont_L = 1.0;
  real cont_H = 0.09;
  real cont_W = 0.09;
  
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
    // return false && (((0.5 - bmarg) < r[0] && r[0] < (0.5 + 0.00625 / 1.0 + bmarg)) &&
    // 	    ((0.0 - bmarg) < r[1] && r[1] < (0.2 + bmarg)) &&
    // 	    ((0.45 - bmarg) < r[2] && r[2] < (0.55 + bmarg)));
    return true;
  }
}

real xmin = 0.0;
real xmax = Geo::xmax;
real ymin = 0.0;
real ymax = Geo::ymax;
real zmin = 0.0;
real zmax = Geo::zmax;
real xcenter = 0.2;
real ycenter = 0.2;
real radius = 0.05;

real density(Point p)
{
  real val = 0.0;

  //val = 1.1 - 0.1 * p[1];
  val = 1.0;

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

// Sub domain for outflow
class OutflowBoundary : public SubDomain
{
public:
  bool inside(const real* p, bool on_boundary) const
  {
    return on_boundary && (p[0] > (xmax - bmarg));
  }
};

// Sub domain for inflow
class InflowBoundary : public SubDomain
{
public:
  bool inside(const real* p, bool on_boundary) const
  {
    return on_boundary && (p[0] < (xmin + bmarg));
  }
};

// Sub domain for inflow
class FixBoundary : public SubDomain
{
public:
  bool inside(const real* r, bool on_boundary) const
  {
    return on_boundary && false &&
      ((r[0] > 0.32 - bmarg) &&
       (r[0] < 0.33 + bmarg) &&
       (r[1] > 0.045 - bmarg) &&
       (r[1] < 2.0 + bmarg) &&
       (r[2] > 0.03- bmarg) &&
       (r[2] < 0.06 + bmarg));
  }
};

// Sub domain for outflow
class SlipBoundary : public SubDomain
{
public:
  bool inside(const real* p, bool on_boundary) const
  {
    InflowBoundary ib;
    OutflowBoundary ob;
    FixBoundary fb;

    // return on_boundary && !ib.inside(p, on_boundary) &&
    //   !ob.inside(p, on_boundary) && !fb.inside(p, on_boundary);
    return on_boundary && !ib.inside(p, on_boundary) &&
      !ob.inside(p, on_boundary);
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
    real tramp = 1.0e-2;

    

    if(t < tramp)
    {
      ramp = (t - 0.0) / (tramp - 0.0);
    }
    else if(t >= tramp)
    {
      ramp = 1.0;
    }

    //ramp = 1.0;


    for(int i = 0; i < 3; i++)
    {
      if (i==0)
      {
	if ( x[0] < (0.0 + bmarg))
	  {
	    real ubar = dolfin_get("Ubar");

	    values[i] = ramp * ubar;
	  }
	else
	{
	  values[i] = 0.0;
	}
      } 
      else
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


void transform(Mesh& mesh)
{
  // Define the area that the pin occupies
  for(VertexIterator v(mesh); !v.end(); ++v)
  {
    Vertex& vertex = *v;
    
    Point r = vertex.point();

    MeshGeometry& geometry = mesh.geometry();

    geometry.x(vertex.index(), 0) *= 1.0e-3;
    geometry.x(vertex.index(), 1) *= 1.0e-3;
    geometry.x(vertex.index(), 2) *= 1.0e-3;
  }
}

void solve(Mesh& mesh, Checkpoint& chkp, long& w_limit, timeval& s_time, Mesh* structure_mesh)
{
  real T = dolfin_get("T");
  real dual_T = dolfin_get("dual_T");
  real nu = dolfin_get("nu");
  real ubar = dolfin_get("Ubar");
  real rho_s = dolfin_get("rho_s");
  real mu = dolfin_get("mu");
  real nu_s = dolfin_get("nu_s");

  TimeDependent td;

  NodeNormal node_normal(mesh);

  InflowBoundary iboundary;
  BC_Momentum_3D bcf_mom(mesh, td);
  BCDensity density_function(mesh);
  Dual_BC_Momentum_3D dual_bcf_mom(mesh);
  BC_Continuity bcf_con(mesh);
  ForceFunction f(mesh);
  Function zero(mesh, 2, 0.0);
  OutflowBoundary oboundary;
  SlipBoundary slipboundary;
  FixBoundary fixboundary;

  Phi phi(mesh);
  Beta beta(mesh);

  Array<Function*> aero_f;
  aero_f.push_back(&phi);

  //NodeNormal node_normal(mesh);
  DirichletBC p_bcin(bcf_mom, mesh, iboundary);
  DirichletBC p_bcfix(bcf_mom, mesh, fixboundary);
  DirichletBC p_bcout(bcf_con, mesh, oboundary);
  DirichletBC p_bcdensity(density_function, mesh, iboundary);
  SlipBC slip_bc(mesh, slipboundary, node_normal);
  //DirichletBC slip_bc(bcf_mom, mesh, slipboundary);
  
  //dolfin_set("PDE slip alpha max", DOLFIN_PI / 1.9);
  Array <BoundaryCondition*> p_bc_momentum;
  p_bc_momentum.push_back(&slip_bc);
  p_bc_momentum.push_back(&p_bcin);
  p_bc_momentum.push_back(&p_bcfix);

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
  
  int svcounter = 0;
  for (VertexIterator v(mesh); !v.end(); ++v)
  {
    Vertex& vertex = *v;
    Point pref = vertex.point();
    
    if (existRegionFile)
    {
      Array<unsigned int> overlap_cells;
      bool bfnd = false;

      for(int ii = -1; ii <= 1; ii++)
      {
	for(int jj = -1; jj <= 1; jj++)
	{
	  for(int kk = -1; kk <= 1; kk++)
	  {
	    Point p(pref[0] + gmarg*ii, pref[1] + gmarg*jj, pref[2] + gmarg*kk);
		    
	    bfnd = false || bfnd;
	    overlap_cells.clear();
	    idetector->overlap(p, overlap_cells);      

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
	}
      }

      if(bfnd && !mesh.distdata().is_ghost(v->index(), 0)) // 6606 ref
	svcounter++;
      
    }
    else
      solid_vertices.set(vertex, Geo::isStructure(pref));
  }
  
  int svcounter_tmp = svcounter;
  MPI_Allreduce(&svcounter_tmp, &svcounter, 1, MPI_INT, 
		MPI_SUM, dolfin::MPI::DOLFIN_COMM);
  std::cout << "svcounter: " << svcounter << std::endl;

  Function U, U0;

//   dolfin_set("Krylov relative tolerance", 1.0e-12);
//   dolfin_set("Krylov absolute tolerance", 1.0e-20);

  NSESolver psolver(mesh, U, U0, f, f, phi, beta, p_bc_momentum, p_bc_pressure,
		    p_bc_density, &density, solid_cells, solid_vertices, T, nu, mu, nu_s, rho_s,
                    ubar, td, "primal"); 
  dolfin_set("Adaptive refinement percentage", 5.0);
  dolfin_set("ODE discrete tolerance", 1.0e-3);
  dolfin_set("ODE maximum iterations", 30);
  dolfin_set("PDE number of samples", 100);  

  psolver.solve(U, U0);
  
}


int main(int argc, char* argv[])
{
  // UnitCube mesh(80, 20, 20);

  // transform(mesh);

  // File meshfile("cube.xml");

  // meshfile << mesh;

  timeval s_time;
  gettimeofday(&s_time, NULL);
  Mesh mesh;  
  long w_limit = 0;
  Checkpoint chkp;
  int iter = 0;
  Mesh* structure_mesh;

  unicorn_init(argc, argv, mesh, chkp, w_limit, iter, structure_mesh);

  transform(mesh);

  mesh.refine();
  //mesh.refine();

  for(int i = 0; i < 2; i++)
  {
    MeshFunction<bool> cell_refinement_marker(mesh);
    cell_refinement_marker.init(mesh.topology().dim());
    
    for (CellIterator c(mesh); !c.end(); ++c)
    {
      Point tp(0.2, 0.045, 0.045);
      Point r = c->midpoint();

      if(((0.2 - bmarg) < r[0] && r[0] < (0.5 + bmarg)) &&
	 ((0.0 - bmarg) < r[1] && r[1] < (0.09 + bmarg)) &&
	 ((0.0 - bmarg) < r[2] && r[2] < (0.09 + bmarg)) &&
	 c->diameter() > 0.004)
      {
	cell_refinement_marker.set(c->index(), true);
      }
      else
      {
	cell_refinement_marker.set(c->index(), false);
      }
    }
    
    if(dolfin::MPI::processNumber() == 0)
      dolfin_set("output destination","terminal");
    
    const std::string refine_type = dolfin_get("adapt_algorithm");
    if(refine_type == "rivara")
      RivaraRefinement::refine(mesh, cell_refinement_marker);
    else if(refine_type == "simple")
      mesh.refine(cell_refinement_marker, true);
    else
      dolfin::error("Unknown refinement algorithm");
    dolfin_set("output destination","silent");
    
    if(dolfin::MPI::processNumber() == 0)
      dolfin_set("output destination","terminal");
    message("cells after: %d", mesh.distdata().global_numCells());
    message("vertices after: %d", mesh.distdata().global_numVertices());
    dolfin_set("output destination","silent"); 
  }
  
  //mesh.refine();

  for(int i = 0; i < 0; i++)
  {
    MeshFunction<bool> cell_refinement_marker(mesh);
    cell_refinement_marker.init(mesh.topology().dim());
    
    for (CellIterator c(mesh); !c.end(); ++c)
    {
      cell_refinement_marker.set(c->index(), true);
    }

    if(dolfin::MPI::processNumber() == 0)
      dolfin_set("output destination","terminal");
    
    const std::string refine_type = dolfin_get("adapt_algorithm");
    if(refine_type == "rivara")
      RivaraRefinement::refine(mesh, cell_refinement_marker);
    else if(refine_type == "simple")
      mesh.refine(cell_refinement_marker, true);
    else
      dolfin::error("Unknown refinement algorithm");
    dolfin_set("output destination","silent");
    
    if(dolfin::MPI::processNumber() == 0)
      dolfin_set("output destination","terminal");
    message("cells after: %d", mesh.distdata().global_numCells());
    message("vertices after: %d", mesh.distdata().global_numVertices());
    dolfin_set("output destination","silent"); 
  }

  transform(*structure_mesh);
  unicorn_solve(mesh, chkp, w_limit, s_time, iter, 0, 0, &solve, structure_mesh);

  dolfin_finalize();
   return 0;
}
