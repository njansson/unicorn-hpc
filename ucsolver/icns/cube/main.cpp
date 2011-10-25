// Copyright (C) 2005 Johan Hoffman.
// Licensed under the GNU GPL Version 2.
//
// Modified by Anders Logg 2005.
// Modified by Niclas Jansson 2008-2010.
// Modified by Johan Jansson 2011
//
//
// Unicorn test problem from "Automated Scientific Computing" of turbulent flow past a cube.
// The output quantity in "aero_f.dat" is the drag force F (second column). The drag coefficient is Cd = 2F/0.2^2.
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

real xmin = -10.0;
real xmax = 30.0;
real ymin = -10.0;
real ymax = 10.0;
real zmin = -10.0;
real zmax = 10.0;
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
    if ( (x[0] > -0.5 - bmarg) && (x[0] < 0.5 + bmarg) &&
	 (x[1] > -0.5 - bmarg) && (x[1] < 0.5 + bmarg) &&
	 (x[2] > -0.5 - bmarg) && (x[2] < 0.5 + bmarg) )
    {
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

    if ( (x[0] > -0.5 - bmarg) && (x[0] < 0.5 + bmarg) &&
	 (x[1] > -0.5 - bmarg) && (x[1] < 0.5 + bmarg) &&
	 (x[2] > -0.5 - bmarg) && (x[2] < 0.5 + bmarg) )
      values[0] = 1.0;
  }
};

class Psi : public Function
{
public:
  Psi(Mesh& mesh) : Function(mesh) {}
  void eval(real* values, const real* x) const
  {
    int d = cell().dim();

    for(int i = 0; i < d; i++)
    {
      values[i] = 0.0;
    }

    if ( (x[0] > -0.5 - bmarg) && (x[0] < 0.5 + bmarg) &&
	 (x[1] > -0.5 - bmarg) && (x[1] < 0.5 + bmarg) &&
	 (x[2] > -0.5 - bmarg) && (x[2] < 0.5 + bmarg) )
      values[1] = 1.0;
  }
};

// Friction coefficient
class Beta : public Function
{
public:
  Beta(Mesh& mesh) : Function(mesh) {}
  void eval(real* values, const real* x) const
  {
    if ( (x[0] > -0.5 - bmarg) && (x[0] < 0.5 + bmarg) &&
	 (x[1] > -0.5 - bmarg) && (x[1] < 0.5 + bmarg) &&
	 (x[2] > -0.5 - bmarg) && (x[2] < 0.5 + bmarg) )
    {
      values[0] = dolfin_get("beta");
    }
    else
    {
      values[0] = 0.0;
    }
  }
};

void solve(Mesh& mesh, Checkpoint& chkp, long& w_limit, timeval& s_time, Mesh* structure_mesh)
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

  Phi phi(mesh);
  Psi psi(mesh);
  Beta beta(mesh);

  Array<Function*> aero_f;
  aero_f.push_back(&phi);

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
 
  Array <BoundaryCondition*> dual_bc_mom;
  dual_bc_mom.push_back(&dual_bc_mom0);
  dual_bc_mom.push_back(&dual_bc_mom1);  
  dual_bc_mom.push_back(&dual_bc_mom2);

  TimeDependent td;

  timeval st;  
  if (w_limit) {
    gettimeofday(&st, NULL);
    w_limit -= (st.tv_sec - s_time.tv_sec);
  }

  NSESolver psolver(mesh, node_normal, f, beta, aero_f, bc_mom, bc_con,
		    T, nu, ubar, chkp, w_limit, td, "primal"); 
  psolver.solve();

  if (w_limit) {
    gettimeofday(&st, NULL);
    w_limit -= (st.tv_sec - s_time.tv_sec);
  }

  NSESolver dsolver(mesh, node_normal, f, beta, aero_f, dual_bc_mom,
  		    dual_bc_con, dual_T, nu, ubar, chkp, w_limit, td, "dual");
  dsolver.solve();
    
  if (w_limit) {
    gettimeofday(&st, NULL);
    w_limit -= (st.tv_sec - s_time.tv_sec);
  }
  
}

int main(int argc, char* argv[])
{
  timeval s_time;
  gettimeofday(&s_time, NULL);
  Mesh mesh;  
  Mesh* structure_mesh = 0;
  long w_limit = 0;
  Checkpoint chkp;
  int iter = 0;

  unicorn_init(argc, argv, mesh, chkp, w_limit, iter, structure_mesh);

  mesh.refine();

  unicorn_solve(mesh, chkp, w_limit, s_time, iter, 0, 0, &solve, structure_mesh);

  dolfin_finalize();
  return 0;
}
