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
#include <unicorn/SlipBC.h>
#include <unicorn/UniParameters.h>
#include <unicorn/NSESolver.h>
 
using namespace dolfin;
using namespace dolfin::unicorn;

real bmarg = 1.0e-3 + DOLFIN_EPS;

real rd = 5.0;
real b = 0.0;
real coef = 0.0;
real T0 = 0.0;;


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

// Sub domain for outflow
class OutflowBoundary : public SubDomain
{
public:
  bool inside(const real* p, bool on_boundary) const
  {
    return false;
  }
};

// Sub domain for inflow
class InflowBoundary : public SubDomain
{
public:
  bool inside(const real* p, bool on_boundary) const
  {
    return on_boundary 
      && (sqr(p[0]) + sqr(p[1]) + sqr(p[2]) > sqr(rd) - bmarg);
  }
};

// Sub domain for Slip boundary condition
class SlipBoundary : public SubDomain
{
public:
  bool inside(const real* p, bool on_boundary) const
  {
    //    return on_boundary && (p[0] > xmin + bmarg && p[0] < xmax - bmarg); 
    return on_boundary && (sqr(p[0]) + sqr(p[1]) + sqr(p[2]) < sqr(rd) - bmarg);
  }
};

// Sub domain for Box
class BBoundary : public SubDomain
{
public:
  bool inside(const real* p, bool on_boundary) const
  {
    //     return on_boundary &&
    return on_boundary && (sqr(p[0]) + sqr(p[1]) + sqr(p[2]) > sqr(rd) - bmarg);
  }
};

class Phi : public Function
{
public:
  Phi(Mesh& mesh, TimeDependent& td) : Function(mesh), td(td) {}
  void eval(real* values, const real* x) const
  {
    real t = td.time();
    
    real aoa = 0.0;

    // the equation of line:
    if (t < T0)
      aoa = 0.0;
    else
      aoa = coef*t + b;

    int d = cell().dim();
    if ( (sqr(x[0]) + sqr(x[1]) + sqr(x[2]) < rd - bmarg) &&
         fabs(x[2]) < 0.5 - bmarg)
      {
        //orthogonal to the drag direction
        values[0] =  cos(M_PI / 180.0 * aoa);
        values[1] =  sin(M_PI / 180.0 * aoa);
        values[2] =  0.0;
      }
    else
      for(int i = 0; i < d; i++)
        values[i] = 0.0;
  }
  TimeDependent& td;
};

class Psi : public Function
{
public:
  Psi(Mesh& mesh, TimeDependent& td) : Function(mesh), td(td) {}
  void eval(real* values, const real* x) const
  {
    real t = td.time();
    
    real aoa = 0.0;

    // the equation of line:
    if (t < T0)
      aoa = 0.0;
    else
      aoa = coef*t + b;

    int d = cell().dim();
    if ( (sqr(x[0]) + sqr(x[1]) + sqr(x[2]) < rd - bmarg) &&
         fabs(x[2]) < 0.5 - bmarg)
    {
      //orthogonal to the drag direction
      values[0] =  sin(M_PI / 180.0 * aoa);
      values[2] = -cos(M_PI / 180.0 * aoa);
      values[1] =  0.0;
    }
    else
      for(int i = 0; i < d; i++)
        values[i] = 0.0;
  }
  TimeDependent& td;
};

// Friction coefficient
class Beta : public Function
{
public:
  Beta(Mesh& mesh) : Function(mesh) {}
  void eval(real* values, const real* x) const
  {
    if ( (sqr(x[0]) + sqr(x[1]) + sqr(x[2]) < rd - bmarg) &&
         fabs(x[2]) < 0.5 - bmarg)
      {
        values[0] = dolfin_get("beta");
      }
    else
      {
        values[0] = 0.0;
      }
  }
};

class TimeDependentBC : public Function
{
public:
  TimeDependentBC(Mesh& mesh, TimeDependent& td) : Function(mesh), td(td) {}
  void eval(real* values, const real* x) const
  {
    real t = td.time();
   
    real aoa = 0.0;

    // the equation of line:
    if (t < T0)
      aoa = 0.0;
    else
      aoa = coef*t + b;
    
    values[0] = cos(M_PI / 180.0 * aoa);
    values[1] = sin(M_PI / 180.0 * aoa);
    values[2] = 0.0;
  }
  TimeDependent& td;

};


class TimeDependentDualBC : public Function
{
public:
  TimeDependentDualBC(Mesh& mesh, TimeDependent& td) : Function(mesh), td(td) {}
  void eval(real* values, const real* x) const
  {
    real t = td.time();
    
    real aoa = 0.0;

    // the equation of line:
    if (t < T0)
      aoa = 0.0;
    else
      aoa = coef*t + b;

    values[0] = - cos(M_PI / 180.0 * aoa);
    values[1] = - sin(M_PI / 180.0 * aoa);
    values[2] = 0.0;
  }
  TimeDependent& td;

};

void solve(Mesh& mesh, Checkpoint& chkp, long& w_limit, timeval& s_time)
{
  
  real T = dolfin_get("T");
  real dual_T = dolfin_get("dual_T");
  real nu = dolfin_get("nu");
  real ubar = dolfin_get("Ubar");

  TimeDependent td;
  Phi phi(mesh, td);
  Psi psi(mesh, td);
  Beta beta(mesh);
  ForceFunction f(mesh);

  Array <Function*> aero_f;
  aero_f.push_back(&phi);
  aero_f.push_back(&psi);
  
  Function zero(mesh, 0.0);
  Function one(mesh, 1.0);
  Function minus_one(mesh, -1.0);
  Function inflow(mesh, ubar);

  InflowBoundary ib;
  OutflowBoundary ob;
  SlipBoundary sb;
  BBoundary bb;

  TimeDependentBC bc_u(mesh, td);
  DirichletBC p_bcin(bc_u, mesh, ib);
  DirichletBC p_bcout(zero, mesh, ob);
  NodeNormal node_normal(mesh);
  SlipBC slip_bc(mesh, sb, node_normal);

  Array <BoundaryCondition*> p_bc_mom;
  p_bc_mom.push_back(&p_bcin);
  p_bc_mom.push_back(&slip_bc);

  TimeDependentDualBC bc_du(mesh, td);
  DirichletBC d_bc1(zero, mesh, bb);
  DirichletBC d_cbc(bc_du, mesh, sb);
  DirichletBC d_bcout(zero, mesh, ob);
  
  Array <BoundaryCondition*> d_bc_mom;
  d_bc_mom.push_back(&d_bc1);
  d_bc_mom.push_back(&d_cbc);

  timeval st;  
  if (w_limit) {
    gettimeofday(&st, NULL);
    w_limit -= (st.tv_sec - s_time.tv_sec);
  }

  NSESolver psolver(mesh, node_normal, f, beta, aero_f, p_bc_mom, p_bcout, 
		    T, nu, ubar, chkp, w_limit, td, "primal"); 
  psolver.solve();

  if (w_limit) {
    gettimeofday(&st, NULL);
    w_limit -= (st.tv_sec - s_time.tv_sec);
  }

  NSESolver dsolver(mesh, node_normal, f, beta, aero_f, d_bc_mom, d_bcout, 
		    dual_T, nu, ubar, chkp, w_limit, td, "dual");
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
  long w_limit = 0;
  Checkpoint chkp;
  int iter = 0;
  unicorn_init(argc, argv, mesh, chkp, w_limit, iter);

  real T = dolfin_get("T");
  real alpha = dolfin_get("alpha");
  T0 = T/10.0;
  coef = alpha/ (T - T0);
  b = - T0 * coef;
  
  unicorn_solve(mesh, chkp, w_limit, s_time, iter, 0, 0, &solve);

  dolfin_finalize();
  return 0;
}
