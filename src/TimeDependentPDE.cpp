// Copyright (C) 2006 Johan Jansson.
// Licensed under the GNU GPL Version 2.
//
// First added:  2006
// Last changed: 2006-05-04


#include <dolfin.h>

#include "unicorn/Project.h"
#include "unicorn/TimeDependentPDE.h"

using namespace dolfin;
using namespace dolfin::unicorn;

//-----------------------------------------------------------------------------
TimeDependentPDE::TimeDependentPDE(Form& a, Form& L, Mesh& mesh, 
				   Array <BoundaryCondition*>& bcs, real T,
				   std::string name,
				   Form* aJac, Form* aP, Form* LP)
  : _mesh(&mesh), _bcs(&bcs),
    t(0), T(T),
    k(dolfin_get("ODE initial time step")),
    solutionfile(dolfin_get("solution file name"), t),
    x(0), dotx(0),
    fevals(0),
    itercounter(0),
    sampleperiod(T / 100.0), lastsample(0.0),
    ksolver(0),
    assembler(0),
    reset_tensor(true),
    reassemble(true),
    J(0), JNoBC(0),
    name(name)
{

  if(!ParameterSystem::parameters.defined("PDE reassemble matrix"))
    dolfin_add("PDE reassemble matrix", reassemble);

  init(a, L, aJac, aP, LP);

}
//-----------------------------------------------------------------------------
TimeDependentPDE::TimeDependentPDE(Mesh& mesh,
				   Array <BoundaryCondition*>& bcs, real T, std::string name) :
    _mesh(&mesh), _bcs(&bcs),
    t(0), T(T),
    k(dolfin_get("ODE initial time step")),
    solutionfile(dolfin_get("solution file name")),
    x(0), dotx(0),
    fevals(0), 
    sampleperiod(T / 100.0), lastsample(0.0),
    assembler(0),
    reset_tensor(true),
    reassemble(true),
    J(0), JNoBC(0),
    name(name)
{

  if(!ParameterSystem::parameters.defined("PDE reassemble matrix"))
    dolfin_add("PDE reassemble matrix", reassemble);

  // Up to subclass to initialize forms
}
//-----------------------------------------------------------------------------
void TimeDependentPDE::init(Form& a, Form& L, Form* aJac,
			    Form* aP, Form* LP)
{
  _a = &a;
  _Lf = &L;
  this->aJac = aJac;
  this->aP = aP;
  this->LP = LP;
}
//-----------------------------------------------------------------------------
TimeDependentPDE::~TimeDependentPDE()
{
  free();
}
//-----------------------------------------------------------------------------
dolfin::uint TimeDependentPDE::solve(Function& U, Function& U0)
{
  dolfin_assert(_a);
  dolfin_assert(_Lf);
  dolfin_assert(_mesh);

  pde_timer = time(); //.start();

  message("Solving time dependent PDE.");

  u0(U.vector());
  U.vector().zero();

  save(U, t);

  // Start time-stepping
  while(t < T) {

    total_timer = time(); //.start();
    
    local_timer = time(); //.start();
    if(assembler)
    {
      delete assembler;
      assembler = new Assembler(mesh());
    }
    else
    {
      assembler = new Assembler(mesh());
    }

    preparestep();
    //message("TPDE preparestep timer: %g", local_timer.elapsed());
    message("TPDE preparestep timer: %g", time() - local_timer);
//    local_timer.stop();

    //local_timer.restart();
    local_timer = time(); //.start();
    *x0 = *x;

    step();
    //message("TPDE step timer: %g", local_timer.elapsed());
    message("TPDE step timer:  %g", time() - local_timer);
//    local_timer.stop();

    local_timer = time(); //.start();
    save(U, t);

    local_timer = time(); //.start();

    bool cont = update(t, false); 
    if(!cont)
      break;

    t += k;
    cout << "TPDE: Stepping k: " << k << " t: " << t << endl;
    cout << "TPDE total step timer: " << time() -  total_timer<< endl;;
  }

  message("TPDE total PDE timer:  %g", time() - pde_timer);
  //pde_timer.stop();

  return 0;
}
//-----------------------------------------------------------------------------
void TimeDependentPDE::u0(GenericVector& x0)
{
  x0.zero();
}
//-----------------------------------------------------------------------------
void TimeDependentPDE::step()
{
  const real tol = dolfin_get("ODE discrete tolerance");

  cout << "TimeDependentPDE::step" << endl;

  incr = 0;
  real incr0 = 0;
  real res0 = 0;
  int maxit = dolfin_get("ODE maximum iterations");
  for(int it = 0; it < maxit; it++)
  {
    incr0 = incr;
    res0 = res;
    res = 0;
    local2_timer = time(); //.start();
    prepareiteration();
    cout << "TPDE prepareiteration timer: " << time() - local2_timer << endl;
    cout << "maxit: " << maxit << endl;
    incr = iter();
    itercounter++;
    postiteration();

    cout << name << " increment: " << incr << " " << it << endl;
    cout << "res: " << res << endl;
    if(incr < tol)
    {
      cout << "increment: iteration converged" << endl;
      break;
    }
    else if(it > maxit)
    {
      message("TPDE: fixed-point iteration diverging");
      k = 0.5 * k;
      reassemble = true;
      *x = *x0;
      revert();
      it = 0;
    }
    else if(it == maxit - 1)
    {
      message("TPDE: fixed-point iteration didn't converge");
    }
    else
    {
      cout << "continue iter" << endl;
    }
    reset_tensor = false;
  }
}
//-----------------------------------------------------------------------------
real TimeDependentPDE::iter()
{
  cout << "TimeDependentPDE::iter" << endl;
  cout << "iterk: " << k << endl;

  iter_timer = time(); //.start();

  local_timer= time(); //.start();
  //if(t == 0.0)

  reassemble = true;

  if(reassemble)
    dolfin_set("PDE reassemble matrix", true);
  else
    dolfin_set("PDE reassemble matrix", false);


  message("TPDE: reassembling Jacobian matrix");
  assembler->assemble(*J, a(), reset_tensor);
  cout << "TPDE assemble J timer: " << time() - local_timer << endl;

  local_timer = time() ;
  fu(*x, *dotx, t);
  cout << "TPDE fu timer: " << time() - local_timer << endl;
//  local_timer.stop();
  // Add J * UP to Newton residual
//   JNoBC->mult(*x, *xtmp);
//   *dotx += *xtmp;


  local_timer = time(); //.start();
  for (uint i = 0; i < bc().size(); i++)
    bc()[i]->apply(*J, *dotx, a());
  cout << "TPDE BC timer: " << time() - local_timer << endl;

  if(reset_tensor)
  {
    reset_tensor = false;
  }

  *dx = *x;

  J->mult(*x, *residual);
  *residual -= *dotx;

  local_timer = time(); //. start();
  ksolver->solve(*J, *x, *dotx);
  //lusolver->solve(J, *dx, *dotx);

//   cout << "J:" << endl;
//   J.disp();
//   cout << "x:" << endl;
//   x->disp();
//   cout << "dotx:" << endl;
//   dotx->disp();

  cout << "TPDE linear solve timer: "<< time() - local_timer<< endl; 

  *dx -= *x;
  
  real relincr = 0.0;
  if(x->norm(l2) > 1.0e-8)
    relincr = dx->norm(l2) / x->norm(l2);
  real resnorm = residual->norm(linf);

  real oldres = res;
  res = std::max(oldres, resnorm);

  cout << "TPDE total iter timer: " << time() - iter_timer << endl;

  return relincr;
  //return dx->norm(linf);
  //return dx->norm(linf);
}
//-----------------------------------------------------------------------------
bool TimeDependentPDE::update(real t, bool end)
{
  return true;
}
//-----------------------------------------------------------------------------
void TimeDependentPDE::revert()
{
}
//-----------------------------------------------------------------------------
void TimeDependentPDE::reset(real T)
{
  this->T = T;
  this->t = 0;
  itercounter = 0;
  sampleperiod = T / 100.0;
  lastsample = 0.0;
}    
//-----------------------------------------------------------------------------
real TimeDependentPDE::timestep(real t, real k0) const
{
  return k;
}
//-----------------------------------------------------------------------------
void TimeDependentPDE::fu(const Vector& foox, Vector& foodotx, real t)
{
  assembler->assemble(*dotx, L(), reset_tensor);
  cout << "DEBUG dotx: " << dotx->norm(linf) << endl;
  cout << "DEBUG U: " << foox.norm(linf) << endl;
  cout << "DEBUG U2: " << x->norm(linf) << endl;

  fevals = fevals + 1;
}
//-----------------------------------------------------------------------------
Form& TimeDependentPDE::a()
{
  dolfin_assert(_a);
  return *_a;
}
//-----------------------------------------------------------------------------
Form& TimeDependentPDE::L()
{
  dolfin_assert(_Lf);
  return *_Lf;
}
//-----------------------------------------------------------------------------
Mesh& TimeDependentPDE::mesh()
{
  dolfin_assert(_mesh);
  return *_mesh;
}
//-----------------------------------------------------------------------------
Array <BoundaryCondition*>& TimeDependentPDE::bc()
{
  dolfin_assert(_bcs);
  return *_bcs;
}
//-----------------------------------------------------------------------------
void TimeDependentPDE::init(Function& U, Function& U0)
{
  message("TPDE::init(U, U0)");

  if(!J)
    J = new Matrix;

  x = new Vector;
  x0 = new Vector;
  U.init(mesh(), *x, a(), 0);
  U0.init(mesh(), *x0, a(), 0);
  xp = new Vector(x->local_size());
  dotx = new Vector(x->local_size());
  xtmp = new Vector(x->local_size());
  dx = new Vector(x->local_size());
  residual = new Vector(x->local_size());

  assembler = new Assembler(mesh());

  ksolver = new KrylovSolver(bicgstab, sor);

//  lusolver = new LUSolver;

}
//-----------------------------------------------------------------------------
void TimeDependentPDE::free()
{
  //cout << "TPDE::free()" << endl;

  delete x;
  delete x0;
  delete xp;
  delete dotx;
  delete xtmp;
  delete dx;
  delete residual;
  delete ksolver;

  if(JNoBC != 0)
    delete JNoBC;
//   if(J != 0)
//     delete J;
  if(assembler)
    delete assembler;
}
//-----------------------------------------------------------------------------
void TimeDependentPDE::save(Function& U, real t)
{
  //cout << "Saving" << endl;
  
  if(t == 0.0)
  {
    solutionfile << U;
  }
  
  while(lastsample + sampleperiod < t)
  {
    lastsample = std::min(t, lastsample + sampleperiod);
    solutionfile << U;
  }
}
void TimeDependentPDE::preparestep()
{
}
//-----------------------------------------------------------------------------
void TimeDependentPDE::prepareiteration()
{
}
//-----------------------------------------------------------------------------
void TimeDependentPDE::postiteration()
{
}
//-----------------------------------------------------------------------------
void TimeDependentPDE::computeJ()
{
}
//-----------------------------------------------------------------------------
void TimeDependentPDE::computeM()
{
}
