// Copyright (C) 2005 Johan Hoffman.
// Licensed under the GNU GPL Version 2.
//
// Modified by Anders Logg 2006.
//
// First added:  2005
// Last changed: 2006-05-07

#ifndef __NSE_SOLVER_H
#define __NSE_SOLVER_H

#include <dolfin.h>
#include <unicorn/ErrorEstimate.h>
#include <unicorn/TimeDependentPDE.h>

namespace dolfin {
  namespace unicorn {
  // This is a solver for the time dependent incompressible 
  // Navier-Stokes equations. 

    class TimeStepFunction;
    
    class NSESolver : public TimeDependentPDE
    {
    public:
      
      // Create the Navier-Stokes solver
      NSESolver(Mesh& mesh, Function& U, Function& U0,
		Function& f, Function& phi, Function& beta,
		Array<BoundaryCondition*>& bc_mom,
		BoundaryCondition& bc_con,
		real T, real nu, real ubar, TimeDependent& td,
		std::string solver_type);
      
      ~NSESolver();

      // Solve Navier-Stokes equations
      void solve_old();
      
      // Compute cell diameter
      static void ComputeCellSize(Mesh& mesh, Vector& hvector);
      
      // Get minimum cell diameter
      static void GetMinimumCellSize(Mesh& mesh, real& hmin);
      
      // Compute stabilization 
      void ComputeStabilization(Mesh& mesh, Function& w, real nu,
				real k, Vector& d1vector, Vector& d2vector,
				Vector& d1invvector, Form& form);
      
      // Compute cell meanb
      void ComputeMean(Mesh& mesh, Function& vc, Function& vm,
		       Function& v, Function& v0, Form& form, Form& form2);
      
      // Set initial velocity 
      void SetInitialVelocity(Vector& xvel);
      
      // Compute the volume inverse
      void ComputeVolInv(Mesh& mesh, Vector& vol_inv, Form& form);

      // Compute the time derivative of w
      void ComputeTimeDerivative(Mesh& mesh, Function& w,
				 Function& w0, real k,
				 Function& dtw);

      void computeP();
      void prepareiteration();
      void postiteration();
      void preparestep();
      bool update(real t, bool end);
      void save(Function& U, real t);
      
    protected:

      Function& U;
      Function& U0;
      
      Function& f;
      Function& phi;
      Function& beta;
      Array<BoundaryCondition*>& bc_mom;
      BoundaryCondition& bc_con;

      Function Um;
      Function Uc;


      Function P;   // pressure
      Function P0;   // pressure
      Function delta1, delta2, delta1inv; // stabilization coefficients
      Function* fnu;
      TimeStepFunction* fk;
      Function vol_inv;

      Vector U0x;
      Vector Umx;
      Vector Ucx;
      Vector Px;
      Vector P0x;
      Vector Presidual;
      Vector delta1x;
      Vector delta2x;
      Vector delta1invx;
      Vector vol_invx;
      

      Form* aM;
      Form* LM;      
      Form* aC;
      Form* LC;

      Matrix PM, PMdummy;
      Vector PMD;
      //Vector Px;
      Vector Pb;

      real T;
      real nu;
      real ubar;
      real kbar;

      std::string solver_type;

      ErrorEstimate* errest;
      ErrorEstimate* perrest;

      File pfile;

      TimeDependent& td;
      
      KrylovSolver pressure_solver;
      KSP ksp_pressure;



      real Presnorm;
      real startup;
      real hmin;


      uint *indices, *c_indices;
	  
    };

    class TimeStepFunction : public Function
    {
    public:
    TimeStepFunction(Mesh& mesh) : Function(mesh) {}
      void eval(real* values, const real* p) const
      {
	values[0] = *k;
      }
      
      real* k;
    };
    
    
  }}

#endif
