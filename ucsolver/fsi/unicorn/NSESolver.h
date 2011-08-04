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
#include <unicorn/ElasticSmoother.h>
#include <unicorn/LaplacianSmoother.h>
#include <unicorn/TimeDependentPDE.h>

//#include <boost/timer.hpp>

namespace dolfin {
  namespace unicorn {
  // This is a solver for the time dependent incompressible 
  // Navier-Stokes equations. 

    class PointerConstantFunction;
    
    class NSESolver : public TimeDependentPDE
    {
    public:
      
      // Create the Navier-Stokes solver
      NSESolver(Mesh& mesh, Function& U, Function& U0,
		Function& f, Function& fc,
		Function& phi, Function& beta,
		Array<BoundaryCondition*>& bc_mom, 
		Array<BoundaryCondition*>& bc_con, 
		Array<BoundaryCondition*>& bc_rho, 
		real (*density)(Point p),
		MeshFunction<bool>& solid_cells,
		MeshFunction<bool>& solid_vertices,
		real T, real nu, real mu, real nu_s, real rho_s, real ubar, TimeDependent& td,
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
      void ComputeVolInv(Vector& vol_invx);

      // Compute the time derivative of w
      void ComputeTimeDerivative(Mesh& mesh, Function& w,
				 Function& w0, real k,
				 Function& dtw);

      void computeP();
      void computeS();
      void computeRho();
      void computeX(Function& XX);
      void computeXinc();
      void computeW(bool solid);
      void prepareiteration();
      void postiteration();
      void preparestep();
      bool update(real t, bool end);
      void save(Function& U, real t);

      void smoothMesh();
      void deform(Function& X);
      void deform_fluid(Function& X);
      void deform_solid(Function& X);
      
    protected:

      Function& U;
      Function& U0;
      
      Function& f;
      Function& phi;
      Function& beta;
      Array<BoundaryCondition*>& bc_mom;
      Array<BoundaryCondition*>& bc_con;
      Array<BoundaryCondition*>& bc_rho;

      Function Um;
      Function Uc;


      Function P;   // pressure
      Function P0;   // pressure
      Function delta1, delta2, delta1inv; // stabilization coefficients
      Function* fnu;
      Function* frho_s;
      PointerConstantFunction* fk;
      PointerConstantFunction* muf;
      PointerConstantFunction* bf;
      real mu;
      real bb;
      real rho_s;

      Function vol_inv;

      Function rho;   // density
      Function rho0;   // density

      Function theta;   // phase

      Function X;   // position
      Function X0;   // position
      Function Xinc;   // position increment
      Function motion;   // position
      Function Xtmp;   // position
      Function Xtmp2;   // position

      Function W;   // mesh velocity
      Function W0;   // mesh velocity
      Function Wm;   // mesh velocity

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
      
      Vector rhox;
      Vector rho0x;

      Vector thetax;

      Vector Xx;
      Vector X0x;
      Vector Xincx;

      Vector Xtmpx;
      Vector Xtmp2x;

      Vector motionx;

      Vector Wx;
      Vector W0x;
      Vector Wmx;

      Form* aM;
      Form* LM;      
      Form* aC;
      Form* LC;
      Form* aR;
      Form* LR;
      Form* aS;
      Form* LS;

      Matrix PM, PMdummy;
      Vector PMD;
      //Vector Px;
      Vector Pb;
      Matrix RhoM;
      Vector Rhob;

      Function S;
      Function S0;
      Function Sdot;
      Function Sdot0;

      Vector S0x, Sx, Sdotx, Sdot0x, Stmp, Sresidual;
      Matrix SM;
      Vector Sm;

      real T;
      real nu;
      real ubar;
      real kbar;

      std::string solver_type;

      ErrorEstimate* errest;
      ErrorEstimate* perrest;

      File pfile;
      File rhofile;
      File wfile;
      File thetafile;

      TimeDependent& td;
      
      KrylovSolver pressure_solver;
      KSP ksp_pressure;

      ElasticSmoother* smoother;
      LaplacianSmoother* lsmoother;

      real Presnorm;
      real startup;
      real hmin;
      real mu_bar;

      MeshQuality* mqual;

      uint *indices, *c_indices;

      std::ofstream outpFile;
      bool has_output;
      dolfin::uint v_output;

      MeshFunction<real> h0;
      MeshFunction<bool>& solid_cells;
      MeshFunction<bool>& solid_vertices;

   //   boost::timer timer0;
   //   boost::timer timer1;
  		real timer0;
      real timer1;
      int smooth_counter;

      MeshFunction<real> meshf_theta;
    };

    class PointerConstantFunction : public Function
    {
    public:
    PointerConstantFunction(Mesh& mesh) : Function(mesh) {}
      void eval(real* values, const real* p) const
      {
	values[0] = *k;
      }
      
      real* k;
    };
    
    
  }}

#endif
