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
    JNoBC(0)
{

  if(!ParameterSystem::parameters.defined("PDE reassemble matrix"))
    dolfin_add("PDE reassemble matrix", reassemble);

  init(a, L, aJac, aP, LP);

}
//-----------------------------------------------------------------------------
TimeDependentPDE::TimeDependentPDE(Mesh& mesh,
				   Array <BoundaryCondition*>& bcs, real T) :
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
    JNoBC(0)
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

  message("Solving time dependent PDE.");

  u0(U.vector());
  U.vector().zero();

  save(U, t);

  korig = k;

  // Start time-stepping
  while(t < T) {

    if(k != korig)
    {
      reassemble = true;
      k = korig;
    }

    MPI::startTimer(total_timer);
    MPI::startTimer(local_timer);
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
    if(dolfin::MPI::processNumber() == 0)
      message("TPDE preparestep timer: %g", MPI::stopTimer(local_timer));

    MPI::startTimer(local_timer);
    *x0 = *x;

    korig = k;
    step();
    if(dolfin::MPI::processNumber() == 0)
      message("TPDE step timer: %g", MPI::stopTimer(local_timer));



    save(U, t);

    MPI::startTimer(local_timer);

    bool cont = update(t, false); 
    if(!cont)
      break;

    t += k;
    if(dolfin::MPI::processNumber() == 0)
    {
      std::cout << "TPDE: Stepping k: " << k << " t: " << t << std::endl;
      std::cout << "TPDE total step timer: " << MPI::stopTimer(total_timer) << std::endl;
    }
  }

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

  dolfin_set("output destination", "terminal");
  if(dolfin::MPI::processNumber() == 0)
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
    MPI::startTimer(local2_timer);
    prepareiteration();
    if(dolfin::MPI::processNumber() == 0)
      std::cout << "TPDE prepareiteration timer: " << MPI::stopTimer(local2_timer) << std::endl;
    if(dolfin::MPI::processNumber() == 0)
      std::cout << "maxit: " << maxit << std::endl;
    incr = iter();
    itercounter++;
    postiteration();

    if(dolfin::MPI::processNumber() == 0)
      std::cout << "increment: " << incr << " " << it << std::endl;
    if(dolfin::MPI::processNumber() == 0)
      std::cout << "res: " << res << std::endl;
    if(incr < tol)
    {
      if(dolfin::MPI::processNumber() == 0)
	std::cout << "increment: iteration converged" << std::endl;
      break;
    }
    else if(it > maxit)
    {
      if(dolfin::MPI::processNumber() == 0)
	message("TPDE: fixed-point iteration diverging");
      k = 0.5 * k;
      reassemble = true;
      *x = *x0;
      revert();
      it = 0;
    }
    else if(it == maxit - 1)
    {
      if(dolfin::MPI::processNumber() == 0)
	message("TPDE: fixed-point iteration didn't converge");
    }
    else
    {
      if(dolfin::MPI::processNumber() == 0)
	std::cout << "continue iter" << std::endl;
    }
    reset_tensor = false;
  }
}
//-----------------------------------------------------------------------------
real TimeDependentPDE::iter()
{
  if(dolfin::MPI::processNumber() == 0)
  {
    std::cout << "TimeDependentPDE::iter" << std::endl;
    std::cout << "iterk: " << k << std::endl;
  }

  MPI::startTimer(iter_timer);

  MPI::startTimer(local_timer);
  //if(t == 0.0)

  reassemble = true;

  if(reassemble)
    dolfin_set("PDE reassemble matrix", true);
  else
    dolfin_set("PDE reassemble matrix", false);


  if(dolfin::MPI::processNumber() == 0)
    message("TPDE: reassembling Jacobian matrix");
  dolfin_set("output destination", "silent");
  assembler->assemble(J, a(), reset_tensor);
  dolfin_set("output destination", "terminal");
  if(dolfin::MPI::processNumber() == 0)
    std::cout << "TPDE assemble J timer: " << MPI::stopTimer(local_timer) << std::endl;

  MPI::startTimer(local_timer);
  fu(*x, *dotx, t);
  if(dolfin::MPI::processNumber() == 0)
    std::cout << "TPDE fu timer: " << MPI::stopTimer(local_timer) << std::endl;

  // Add J * UP to Newton residual
//   JNoBC->mult(*x, *xtmp);
//   *dotx += *xtmp;


  MPI::startTimer(local_timer);
  dolfin_set("output destination", "silent");
  for (uint i = 0; i < bc().size(); i++)
    bc()[i]->apply(J, *dotx, a());
  if(dolfin::MPI::processNumber() == 0)
    std::cout << "TPDE BC timer: " << MPI::stopTimer(local_timer) << std::endl;

  if(reset_tensor)
  {
    reset_tensor = false;
  }

  *dx = *x;

  J.mult(*x, *residual);
  *residual -= *dotx;

  MPI::startTimer(local_timer);
  dolfin_set("output destination", "silent");
  if(dolfin::MPI::processNumber() == 0)
      dolfin_set("output destination", "terminal");
  ksolver->solve(J, *x, *dotx);
  dolfin_set("output destination", "silent");
  //lusolver->solve(J, *dx, *dotx);

//   std::cout << "J:" << std::endl;
//   J.disp();
//   std::cout << "x:" << std::endl;
//   x->disp();
//   std::cout << "dotx:" << std::endl;
//   dotx->disp();

  if(dolfin::MPI::processNumber() == 0)
    std::cout << "TPDE linear solve timer: " << MPI::stopTimer(local_timer) << std::endl;

  *dx -= *x;
  
  real relincr = dx->norm(linf) / x->norm(linf);
  real resnorm = residual->norm(linf);

  real oldres = res;
  res = std::max(oldres, resnorm);

  if(dolfin::MPI::processNumber() == 0)
    std::cout << "TPDE total iter timer: " << MPI::stopTimer(iter_timer) << std::endl;

  return relincr;
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
real TimeDependentPDE::timestep(real t, real k0) const
{
  return k;
}
//-----------------------------------------------------------------------------
void TimeDependentPDE::fu(const Vector& foox, Vector& foodotx, real t)
{
  dolfin_set("output destination", "silent");
  assembler->assemble(*dotx, L(), reset_tensor);
  dolfin_set("output destination", "terminal");

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

  ksolver = new KrylovSolver(bicgstab, jacobi);

//  lusolver = new LUSolver;

}
//-----------------------------------------------------------------------------
void TimeDependentPDE::free()
{
  //std::cout << "TPDE::free()" << std::endl;

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
  if(assembler)
    delete assembler;
}
//-----------------------------------------------------------------------------
void TimeDependentPDE::save(Function& U, real t)
{
  //std::cout << "Saving" << std::endl;
  
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
