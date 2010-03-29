// Copyright (C) 2005 Johan Hoffman.
// Licensed under the GNU GPL Version 2.
//
// Modified by Anders Logg 2005.
// Modified by Niclas Jansson 2008-2010.
//
// First added:  2002-11-29
// Last changed: 2010-03-29
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
#include <unicorn/SlipBC.h>
#include <unicorn/UniParameters.h>
#include <unicorn/NSESolver.h>
 
using namespace dolfin;
using namespace dolfin::unicorn;

real bmarg = 1.0e-3 + DOLFIN_EPS;

real xmin = 0.0;
real xmax = 2.1;
real ymin = 0.0;
real ymax = 1.4;
real zmin = 0.0;
real zmax = 0.4;
real xcenter = 0.5;
real ycenter = 0.7;
real radius = 0.05;

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
class InflowTopBottomBoundary3D : public SubDomain
{
public:
  bool inside(const real* p, bool on_boundary) const
  {
    return on_boundary && p[0] < xmin + bmarg ;
  }
};

// Sub domain for slip
class SlipBoundary : public SubDomain
{
public:
  bool inside(const real* p, bool on_boundary) const
  {

    return on_boundary && p[0] < (xmax - bmarg) && p[0] > (xmin + bmarg);
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
  BC_Momentum_3D(Mesh& mesh) : Function(mesh) {}
  void eval(real* values, const real* x) const
  {
    if(x[0] < xmin + bmarg) {
      values[0] = 1.0;
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

void solve(Mesh& mesh, Checkpoint& chkp, long& w_limit, timeval& s_time)
{
  
  real T = dolfin_get("T");
  real dual_T = dolfin_get("dual_T");
  real nu = dolfin_get("nu");
  real ubar = dolfin_get("Ubar");
  
  InflowTopBottomBoundary3D iboundary;
  BC_Momentum_3D bcf_mom(mesh);
  Dual_BC_Momentum_3D dual_bcf_mom(mesh);
  BC_Continuity bcf_con(mesh);
  ForceFunction f(mesh);
  Function zero(mesh, 3, 0.0);
  OutflowBoundary oboundary;
  SlipBoundary slipboundary;
  SideWallBoundary swboundary;

  Phi phi(mesh);
  Beta beta(mesh);

  NodeNormal node_normal(mesh);
  SlipBC slip_bc0(mesh, slipboundary, node_normal);
  DirichletBC bc_mom0(bcf_mom, mesh, iboundary);
  DirichletBC bc_con(bcf_con, mesh, oboundary);

  Array <BoundaryCondition*> bc_mom;
  bc_mom.push_back(&bc_mom0);
  bc_mom.push_back(&slip_bc0);

  DirichletBC dual_bc_con(bcf_con, mesh, oboundary);
  DirichletBC dual_bc_mom0(zero, mesh, iboundary);
  DirichletBC dual_bc_mom1(zero, mesh, oboundary);
  DirichletBC dual_bc_mom2(dual_bcf_mom, mesh, slipboundary);
  DirichletBC dual_bc_mom3(zero, mesh, swboundary); 
 
  Array <BoundaryCondition*> dual_bc_mom;
  dual_bc_mom.push_back(&dual_bc_mom0);
  dual_bc_mom.push_back(&dual_bc_mom1);  
  dual_bc_mom.push_back(&dual_bc_mom2);
  dual_bc_mom.push_back(&dual_bc_mom3);

  TimeDependent td;

  timeval st;  
  if (w_limit) {
    gettimeofday(&st, NULL);
    w_limit -= (st.tv_sec - s_time.tv_sec);
  }

  NSESolver psolver(mesh, node_normal, f, phi, beta, bc_mom, bc_con,
		    T, nu, ubar, chkp, w_limit, td, "primal"); 
  psolver.solve();

  if (w_limit) {
    gettimeofday(&st, NULL);
    w_limit -= (st.tv_sec - s_time.tv_sec);
  }

  NSESolver dsolver(mesh, node_normal, f, phi, beta, dual_bc_mom, dual_bc_con,
		    dual_T, nu, ubar, chkp, w_limit, td, "dual");
  dsolver.solve();
    
  if (w_limit) {
    gettimeofday(&st, NULL);
    w_limit -= (st.tv_sec - s_time.tv_sec);
  }
  
}

void smooth(Mesh& mesh)
{  
  BoundaryMesh b(mesh);
  if(!b.numCells()) return;
  MeshFunction<dolfin::uint>* cell_map = b.data().meshFunction("cell map");  
  for (CellIterator bf(b); !bf.end(); ++bf) 
  {
    Facet f(mesh, cell_map->get(*bf));
    
    for(VertexIterator v(f); !v.end(); ++v)
    {
      real x = v->point().x();
      real y = v->point().y();
           
      if ( (sqr(x - xcenter) + sqr(y - ycenter)) <= sqr(radius + 1e-3) )
      {
        real l = radius - sqrt(sqr(x - xcenter) + sqr(y - ycenter));

        real alpha = atan( - (y - ycenter)/(x - xcenter));
        real dx = l*cos(alpha);
        real dy = l*sin(alpha);

        if (fabs(l) > DOLFIN_EPS)
        {
          dx = fabs(dx)/fabs(l)*l;
          dy = fabs(dy)/fabs(l)*l;
        }

        (x - xcenter < 0.0 ? x -= dx : x += dx);
        (y - ycenter < 0.0 ? y -= dy : y += dy);

        *(mesh.coordinates()+(3 * v->index())) = x;
        *(mesh.coordinates()+(3 * v->index() + 1)) = y;
      }
    }
  }
}

int main(int argc, char* argv[])
{
  timeval s_time;
  gettimeofday(&s_time, NULL);
  Mesh mesh;  
  long w_limit = 0;
  Checkpoint chkp;
  int iter = 0;

  unicorn_init(argc, argv, mesh, chkp, w_limit, iter);

  unicorn_solve(mesh, chkp, w_limit, s_time, iter, 0, &smooth, &solve);

  dolfin_finalize();
  return 0;
}
