// Copyright (C) 2006 Johan Jansson.
// Licensed under the GNU GPL Version 2.
//
// Modified by Garth N. Wells 2006.
//
// First added:  2005
// Last changed: 2006-08-21

#ifndef __TIME_DEPENDENT_PDE_H
#define __TIME_DEPENDENT_PDE_H

#include <dolfin.h>

namespace dolfin { namespace unicorn
{

  /// This class implements the solution functionality for time dependent PDEs.

  class TimeDependentPDE
  {
  public:

    /// Define a time dependent PDE with boundary conditions
    TimeDependentPDE(Form& a, Form& L, Mesh& mesh,
		     Array <BoundaryCondition*>& bcs, real T,
		     Form* aJac = 0,
		     Form* aP = 0,
		     Form* LP = 0);

    TimeDependentPDE(Mesh& mesh, Array <BoundaryCondition*>& bcs, real T);

    /// Destructor
    virtual ~TimeDependentPDE();

    /// Initialize PDE
    virtual void init(Form& a, Form& L, Form* aJac, Form* aP, Form* LP);

    /// Solve PDE (in general a mixed system)
    virtual uint solve(Function& U, Function& U0);

    /// Compute initial value
    virtual void u0(GenericVector& u);

    /// Update PDE
    virtual bool update(real t, bool end);

    /// Revert one time step (e.g. if diverging)
    virtual void revert();

    /// ODE timestep
    virtual real timestep(real t, real k0) const;

    /// Compute right hand side dotu = f(u)
    virtual void fu(const Vector& x, Vector& dotx, real t);

    virtual void init(Function& U, Function& U0);

    virtual void free();

    virtual void save(Function& U, real t);

    void step();
    real iter();

    virtual void preparestep();
    virtual void prepareiteration();
    virtual void postiteration();

    virtual void computeJ();
    virtual void computeM();

    /// Return the bilinear form a(.,.)
    Form& a();

    /// Return the linear form L(.,.)
    Form& L();

    /// Return the mesh
    Mesh& mesh();

    /// Return the boundary condition
    Array <BoundaryCondition*>& bc();

  protected:

    Form* _a;
    Form* _Lf;
    Mesh* _mesh;
    Array <BoundaryCondition*>* _bcs;

    int N;
    real t;
    real T;

    real k;
    real korig;

    File solutionfile;

    real incr;
    real res;

  public:
    Vector* x;
    Vector* xp;
    Vector* x0;
    Vector* dotx;
    Vector* xtmp;
    Vector* dx;
    Vector* residual;

    bool implicit;
    bool constantM;
    bool constantJ;

    int fevals;
    int itercounter;
    real sampleperiod;
    real lastsample;

    KrylovSolver* ksolver;
    LUSolver* lusolver;

    Assembler* assembler;
    Form* aJac;
    Form* aP;
    Form* LP;

    bool reset_tensor;
    bool reassemble;

    Matrix J;
    Matrix* JNoBC;
    Vector JD;

    real local_timer;
    real local2_timer;
    real iter_timer;
    real total_timer;
  };

}}
#endif
