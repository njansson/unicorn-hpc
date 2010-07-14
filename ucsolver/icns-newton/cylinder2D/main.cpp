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
  BC_Momentum_3D(Mesh& mesh) : Function(mesh) {}
  void eval(real* values, const real* x) const
  {
    for(int i = 0; i < 2; i++)
    {
      if (i==0)
      {
	if ( x[0] < (0.0 + bmarg) &&
	     x[1] < (ymax - bmarg) &&
	     x[1] > (ymin + bmarg))
	  {
	    real Um = 1.5;


	    values[i] = 4.0*Um*x[1]*(ymax - x[1]) / (ymax*ymax);
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
  Function zero(mesh, 2, 0.0);
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
  //NewSlipBC slip_bc(mesh, slipboundary, node_normal);
  
  //dolfin_set("PDE slip alpha max", DOLFIN_PI / 1.9);
  Array <BoundaryCondition*> p_bc_mom;
  p_bc_mom.push_back(&p_bcin);
  //p_bc_mom.push_back(&slip_bc);

  TimeDependent td;

  Function U, U0;

//   dolfin_set("Krylov relative tolerance", 1.0e-12);
//   dolfin_set("Krylov absolute tolerance", 1.0e-20);

  NSESolver psolver(mesh, U, U0, f, phi, beta, p_bc_mom, p_bcout, T, nu,
                    ubar, td, "primal"); 
  dolfin_set("Adaptive refinement percentage", 5.0);
  dolfin_set("ODE discrete tolerance", 1.0e-2);
  dolfin_set("PDE number of samples", 100);  

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

  unicorn_init(argc, argv, mesh, chkp, w_limit, iter);

  unicorn_solve(mesh, chkp, w_limit, s_time, iter, 0, 0, &solve);

  dolfin_finalize();
  return 0;
}
