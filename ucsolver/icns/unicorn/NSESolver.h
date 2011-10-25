// Copyright (C) 2005 Johan Hoffman.
// Licensed under the GNU GPL Version 2.
//
// Modified by Anders Logg 2006.
// Modified by Niclas Jansson 2008-2010.
//
// First added:  2005
// Last changed: 2010-10-25

#ifndef __NSE_SOLVER_H
#define __NSE_SOLVER_H

#include <signal.h>
#include <sys/time.h>

#include <dolfin.h>
#include "unicorn/AdaptiveRefinement.h"
#include <unicorn/ErrorEstimate.h>
#include <unicorn/NodeNormal.h>

namespace dolfin {
  namespace unicorn {
  // This is a solver for the time dependent incompressible 
  // Navier-Stokes equations. 
    class NSESolver 
    {
    public:
      
      // Create the Navier-Stokes solver
      NSESolver(Mesh& mesh, NodeNormal& node_normal,
		Function& f, Function& beta,
		Array<Function*>& aero_f,
		Array<BoundaryCondition*>& bc_mom,
		BoundaryCondition& bc_con,
		real T, real nu, real ubar,
		Checkpoint& chkp, long& w_limit,
		TimeDependent& td,
		std::string solver_type);
      
      ~NSESolver();
      
      // Solve Navier-Stokes equations
      void solve();

    private:

      // Cleanup
      void clear();
      
      
      // Get minimum cell diameter
      static void GetMinimumCellSize(Mesh& mesh, real& hmin);
      
      // Compute stabilization 
      void ComputeStabilization(Mesh& mesh, Function& w, real nu,
				       real k, 
				       Vector& d1vector, Vector& d2vector,
				       Form& form);
      
      // Compute cell mean
      void ComputeMean(Mesh& mesh, Function& vmean, Function& v, 
		       Form& form);
      
      // Compute the volume inverse
      void ComputeVolInv(Mesh& mesh, Vector& vol_inv);

      // Compute the time derivative of w
      void ComputeTimeDerivative(Mesh& mesh, Function& w,
				 Function& w0, real k,
				 Function& dtw);

      void ComputeTangentialVectors(Mesh& mesh, 
				    Vector& tau_1,
				    Vector& tau_2,
				    Vector& normal,
				    Form& form,
				    NodeNormal& node_normal);

      // Signal handler
      static void sighandler(int sig_code) 
      {
	switch(sig_code)
	{
	case SIGALRM: 
	  WALL_CLOCK_LIMIT = true; 
	  if(dolfin::MPI::processNumber() == 0)
	    dolfin_set("output destination","terminal");
	  warning("Wall clock limit reached");
	  itimerval itv;
	  itv.it_value.tv_sec = 1800;
	  itv.it_value.tv_usec = itv.it_interval.tv_sec = itv.it_interval.tv_usec = 0;
	  if(setitimer(ITIMER_REAL, &itv, 0) < 0)
	    perror("setitimer failed");
	  dolfin_set("output destination","silent");
	  break;
	default: break;
	}
      }
      
      static SolverType krylov_method(std::string type)
      {
	if (type == "cg")
	  return dolfin::cg;
	else if (type == "gmres")
	  return dolfin::gmres;
	else if (type == "bicgstab")
	  return dolfin::bicgstab;
	else 
	{
	  warning("Undefined solver type            "
		  "Fallback to default krylov method");
	  return dolfin::default_solver;
	}
      }

      static PreconditionerType pc_type(std::string type)
      {
	if (type == "none") 
	  return none;
	else if (type == "bjacobi")
	  return dolfin::default_pc; 
	else if (type == "sor")
	  return dolfin::sor;
	else if (type == "ilu")
	  return dolfin::ilu;
	else if (type == "amg")
	  return dolfin::amg;
	else
	{
	  warning("Undefined preconditioner          "
		  "Fallback to default preconditioner");
	  return default_pc;
	}
      }

      
      Mesh& mesh;
      NodeNormal& node_normal;
      Function& f;
      Function& beta;
      Array<Function*>& aero_f;
      Array<BoundaryCondition*>& bc_mom;
      BoundaryCondition& bc_con;

      real T;
      real nu;
      real ubar;
      real hmin;
      bool schur;

      Checkpoint& chkp;
      long& w_limit;
      TimeDependent& td;
      std::string solver_type;

      ErrorEstimate* errest;
      ErrorEstimate* perrest;

      uint *indices, *c_indices;            
      
      static bool WALL_CLOCK_LIMIT;

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

