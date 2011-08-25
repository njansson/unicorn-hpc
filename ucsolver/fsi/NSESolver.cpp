// Copyrightx (C) 2005 Johan Hoffman.
// Licensed under the GNU GPL Version 2.
//
// Modified by Garth N. Wells 2005.
// Modified by Anders Logg 2005-2006.
//
// First added:  2005
// Last changed: 2006-10-19

#include <sstream>
#include <vector>
#include <algorithm>

#include <dolfin.h>
#include <dolfin/common/timing.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/fem/Assembler.h>
#include <dolfin/la/PETScMatrix.h>
#include <dolfin/la/PETScVector.h>


#include <unicorn/SpaceTimeFunction.h>
#include <unicorn/NodeNormal.h>
#include <unicorn/NSESolver.h>
#include <unicorn/Drag3D.h>

#include <dolfin/la/PETScMatrix.h>
#include <dolfin/la/PETScVector.h>
#include <dolfin/la/PETScKrylovMatrix.h>
#include <dolfin/la/PETScLUSolver.h>

#include <unicorn/NSEMomentum3D.h>
#include "unicorn/NSEContinuity3D.h"
#include <unicorn/NSEDensity3D.h>
#include "unicorn/NSEResMomentum3D.h"
#include "unicorn/NSEResContinuity3D.h"
#include "unicorn/Drag3D.h"

#include "unicorn/NSEDualMomentum3D.h"
#include "unicorn/NSEDualContinuity3D.h"
#include "unicorn/NSEDualGradient3D.h"

#include <unicorn/NSEMomentum2D.h>
#include "unicorn/NSEContinuity2D.h"
#include <unicorn/NSEDensity2D.h>
#include "unicorn/Drag2D.h"
#include <unicorn/NavierStokesStress2D.h>
#include <unicorn/NavierStokesStress3D.h>

#include <ufc.h>
#include <dolfin/fem/UFC.h>
#include <dolfin/fem/Form.h>

//#include <boost/timer.hpp>


using namespace dolfin;
using namespace unicorn;
//-----------------------------------------------------------------------------
NSESolver::NSESolver(Mesh& mesh, Function& U, Function& U0,
		     Function& f, Function& fc,
		     Function& phi, Function& beta,
		     Array<BoundaryCondition*>& bc_mom, 
		     Array<BoundaryCondition*>& bc_con, 
		     Array<BoundaryCondition*>& bc_rho, 
		     real (*density)(Point p),
		     MeshFunction<bool>& solid_cells,
		     MeshFunction<bool>& solid_vertices,
		     real T, real nu, real mu, real nu_s, real rho_s, real ubar, TimeDependent& td,
		     std::string solver_type)
  : TimeDependentPDE(mesh, bc_mom, T, "UCSolver"), U(U), U0(U0),
    f(f), phi(phi), beta(beta),
    bc_mom(bc_mom), bc_rho(bc_rho), bc_con(bc_con),
    T(T), nu(nu), ubar(ubar),
    solver_type(solver_type), errest(0), perrest(0),
    pfile("pressure.pvd"),
    rhofile("density.pvd"),
    wfile("meshvel.pvd"),
    thetafile("theta.pvd"),
    td(td),
    pressure_solver(bicgstab, sor),
    ksp_pressure(0),
    startup(true), indices(0), c_indices(0),
    solid_cells(solid_cells), solid_vertices(solid_vertices),
    smoother(0), lsmoother(0), rho_s(rho_s),
    smooth_counter(0), has_output(false)
{
  // Declare parameters
  
  if(!ParameterSystem::parameters.defined("Adaptive refinement percentage"))
    dolfin_add("Adaptive refinement percentage", 5.0);
  if(!ParameterSystem::parameters.defined("PDE number of samples"))
    dolfin_add("PDE number of samples", 100);

  if(dolfin::MPI::processNumber() == 0)
    dolfin_set("output destination", "terminal");
  else
    dolfin_set("output destination", "silent");

  td.sync(&t);

  // Get the number of space dimensions of the problem 
  int nsd = mesh.topology().dim();

  message("Number of space dimensions: %d",nsd);
  message("Number of global vertices: %d",mesh.distdata().global_numVertices());

  h0.init(mesh, mesh.topology().dim());
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    Cell& cell = *c;

    h0.set(cell, MeshQuality::myDiameter(cell));
  }


  // Set time step (proportional to the minimum cell diameter) 
  GetMinimumCellSize(mesh, hmin);  

  // Take very conservative time-step for startup
  k = 0.3*hmin/ubar;
  message("nu: %f",nu);
  message("ubar: %f",ubar);
  message("hmin: %f",hmin);
  message("k: %f",k);

  fnu = new Function(mesh, nu);
  frho_s = new Function(mesh, rho_s);

  fk = new PointerConstantFunction(mesh);
  fk->k = &k;

  this->mu = mu;
  this->bb = nu_s;

  //mu = 1.0*2.0e9;
  //mu = 2.0e0;

  //bb = sqrt(1.0*mu);

  muf = new PointerConstantFunction(mesh);
  muf->k = &(this->mu);
  Function* lmbdaf = new Function(mesh, 0.0 * 2.0 / 3.0 * this->mu);

  bf = new PointerConstantFunction(mesh);
  bf->k = &bb;

  if ( nsd == 3) 
  {
    if(solver_type == "primal")
    {
      aM = new NSEMomentum3DBilinearForm(Uc, Um, *fnu, delta1, delta2, *fk,
					 rho, theta, *muf, *lmbdaf, W, Wm,
					 *bf, *frho_s);
      LM = new NSEMomentum3DLinearForm(U, U0, Uc, Um, P, 
					 *fnu, delta1, delta2, f, *fk,
				       rho, theta, S, *muf, *lmbdaf, W, Wm, *bf, *frho_s);
      aC = new NSEContinuity3DBilinearForm(delta1, *fk);
      LC = new NSEContinuity3DLinearForm(P0, U, *fk);

      aR = new NSEDensity3DBilinearForm(U, Um, *fk, delta1);
      LR = new NSEDensity3DLinearForm(rho0, U, Um, *fk, delta1);

      aS = new NavierStokesStress3DBilinearForm;
      LS = new NavierStokesStress3DLinearForm(S, U, *muf, vol_inv);
    }
    else if(solver_type == "dual")
    {
      error("Not implemented yet");
    }
  }
  else 
  {
    if(solver_type == "primal")
    {
      aM = new NSEMomentum2DBilinearForm(U, Um, *fnu, delta1, delta2, *fk,
					 rho, theta, *muf, *lmbdaf, W, Wm,
					 *bf);
      LM = new NSEMomentum2DLinearForm(U, U0, Uc, Um, P, 
					 *fnu, delta1, delta2, f, *fk,
				       rho, theta, S, *muf, *lmbdaf, W, Wm);
      aC = new NSEContinuity2DBilinearForm(delta1, *fk);
      LC = new NSEContinuity2DLinearForm(P0, U, *fk);

      aR = new NSEDensity2DBilinearForm(U, Um, *fk, delta1);
      LR = new NSEDensity2DLinearForm(rho0, U, Um, *fk, delta1);

      aS = new NavierStokesStress2DBilinearForm;
      LS = new NavierStokesStress2DLinearForm(S, U, *muf, vol_inv);
    }
    else if(solver_type == "dual")
    {
      error("Not implemented yet");
    }
  }

  Uc.init(mesh, Ucx, *aM, 0);
  Um.init(mesh, Umx, *LM, 4);

  X.init(mesh, Xx, *aM, 0);
  X0.init(mesh, X0x, *aM, 0);
  Xinc.init(mesh, Xincx, *aM, 0);

  motion.init(mesh, motionx, *aM, 0);
  Xtmp.init(mesh, Xtmpx, *aM, 0);
  Xtmp2.init(mesh, Xtmp2x, *aM, 0);

  W.init(mesh, Wx, *aM, 0);
  W0.init(mesh, W0x, *aM, 0);
  Wm.init(mesh, Wmx, *LM, 4);

  P.init(mesh, Px, *aC, 0);
  P0.init(mesh, P0x, *aC, 0);
  delta1.init(mesh, delta1x, *LM, 7);
  delta2.init(mesh, delta2x, *LM, 7);
  delta1inv.init(mesh, delta1invx, *LM, 7);

  rho.init(mesh, rhox, *aR, 0);
  rho0.init(mesh, rho0x, *aR, 0);

  S.init(mesh, Sx, *aS, 0);
  S0.init(mesh, S0x, *aS, 0);
  Sdot.init(mesh, Sdotx, *aS, 0);
  Sdot0.init(mesh, Sdot0x, *aS, 0);
  Stmp.init(S.vector().size());
  Sresidual.init(S.vector().size());

  vol_inv.init(mesh, vol_invx, *LS, 4);

  theta.init(mesh, thetax, *LM, 12);

  Presidual.init(P.vector().local_size());

  P.vector().zero();
  S.vector().zero();
  Sdot.vector().zero();
  W.vector().zero();
  Wm.vector().zero();

  // FIXME: Need to initialize rho
  rho.vector() = 1.0;

  UFC ufc(aR->form(), mesh, aR->dofMaps());
  Cell c(mesh, 0);
  uint local_dim = c.numEntities(0);
  uint *idx  = new uint[nsd * local_dim];
  uint *id  = new uint[nsd * local_dim];
  real *rho_block = new real[nsd * local_dim];  

  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    ufc.update(*cell, mesh.distdata());
    
    (aR->dofMaps())[1].tabulate_dofs(idx, ufc.cell, cell->index());
    
    uint ii = 0;
    uint jj = 0;    
    for(VertexIterator v(*cell); !v.end(); ++v, ii++) 
    {
      Point p = v->point();
      real rho_val = (*density)(p);

      if (!mesh.distdata().is_ghost(v->index(), 0)) 
      {
  	rho_block[jj] = rho_val;
  	id[jj++] = idx[ii];
      }
    }

    //rho.vector().set(rho_block, jj, id);
  }

  //rho.vector().apply();
  delete[] rho_block;
  delete[] idx;
  delete[] id;

  // Initialize X
  {
    int d = mesh.topology().dim();
    UFC ufc(aM->form(), mesh, aM->dofMaps());
    Cell c(mesh, 0);
    uint local_dim = c.numEntities(0);
    uint *idx  = new uint[d * local_dim];
    uint *id  = new uint[d * local_dim];
    real *X_block = new real[d * local_dim];  
    
    for (CellIterator cell(mesh); !cell.end(); ++cell)
    {
      ufc.update(*cell, mesh.distdata());
      (aM->dofMaps())[0].tabulate_dofs(idx, ufc.cell, cell->index());
      
      uint ii = 0;
      uint jj = 0;    
      for(uint i = 0; i < d; i++) 
      {
	for(VertexIterator v(*cell); !v.end(); ++v, ii++) 
	{
	  if (!mesh.distdata().is_ghost(v->index(), 0)) 
	  {
	    X_block[jj] = v->x()[i];
	    id[jj++] = idx[ii];
	  }
	  else
	  {
	  }
	}
      }
      X.vector().set(X_block, jj, id);
    }
    X.vector().apply();
    
    delete[] X_block;
    delete[] idx;
    delete[] id;
  }
  
  meshf_theta.init(mesh, mesh.topology().dim());

  // Initialize theta
  {

    int d = mesh.topology().dim();
    int N = mesh.numVertices();
    if(MPI::numProcesses() > 1)
      N = mesh.distdata().global_numVertices();
    int M = mesh.numCells();
    if(MPI::numProcesses() > 1)
      M = mesh.distdata().global_numCells();
    
    {
      UFC ufc(LM->form(), mesh, LM->dofMaps());
      Cell c(mesh, 0);
      uint local_dim = c.numEntities(0);
      uint *idx  = new uint[local_dim];
      uint *id  = new uint[local_dim];
      real *theta_block = new real[local_dim];  
      real theta_val;

      for (CellIterator cell(mesh); !cell.end(); ++cell)
      {
	ufc.update(*cell, mesh.distdata());
	(LM->dofMaps())[7].tabulate_dofs(idx, ufc.cell, cell->index());
	
	// Only one dof
	uint ii = 0;
	uint jj = 0;    
	if(solid_cells.get(*cell))
	{
	  theta_val = 0.0;
	}
	else
	  theta_val = 1.0;

	theta_block[jj] = theta_val;
	id[jj++] = idx[ii];
	theta.vector().set(theta_block, jj, id);

	meshf_theta.set(*cell, theta_val);
      }
      
      theta.vector().apply();

      delete[] theta_block;
      delete[] idx;
      delete[] id;
    }
  }

  mqual = new MeshQuality(mesh);
  mqual->meshQuality();

  mu_bar = mqual->mu_min;

  // Find output vertex                                                                                               
  for (VertexIterator v(mesh); !v.end(); ++v)
  {
    int vidx = v->index();
    Point p = v->point();

    Point pvout(0.6, 0.21, 0);

    if(p.distance(pvout) < 1.0e-4)
    {
      std::cout << "has_output true" << std::endl;
      v_output = vidx;
      has_output = true;
    }
  }

  // cout << "theta: " << endl;
  std::string ofname = "outp_file.m";
  outpFile.open(ofname.c_str());
  outpFile.close();


  // theta.vector().disp();

  // Initialize vertex mesh function marking solid/fluid vertices

  init(*aM, *LM, 0, 0, 0);

  init(U, U0);
}
//-----------------------------------------------------------------------------
NSESolver::~NSESolver()
{
}
//-----------------------------------------------------------------------------
void NSESolver::save(Function& U, real t)
{
  int nsamples = dolfin_get("PDE number of samples");

   sampleperiod = T / (real)nsamples;

//   //TimeDependentPDE::save(U, t);
   //cout << "Saving" << endl;
  
   std::vector<std::pair<Function*, std::string> > output;
   std::pair<Function*, std::string> U_output(&U, "Velocity");
//    std::pair<Function*, std::string> P_output(&P, "Pressure");
   output.push_back(U_output);
//    output.push_back(P_output);

   if(true || t > 0.1)
   {
     if(t == 0.0)
     {
       solutionfile << output;
       pfile << P;
       rhofile << rho;
       wfile << W;
       thetafile << meshf_theta;
     }
     
     while(lastsample + sampleperiod < t)
     //if(true)
     {
       lastsample = std::min(t, lastsample + sampleperiod);
       solutionfile << output;
       pfile << P;
       rhofile << rho;
       wfile << W;
       thetafile << meshf_theta;
     }
   }


   if(has_output)
   {
     std::string ofname = "outp_file.m";
     outpFile.open(ofname.c_str(), std::ios::app);
     
     Vertex vout(mesh(), v_output);
     real outpval = vout.point()[1];
     
     outpFile << t  << "\t" << outpval << "\n";
     outpFile.flush();
     outpFile.close();
   }
}
//-----------------------------------------------------------------------------
void NSESolver::preparestep()
{
  if(dolfin::MPI::processNumber() == 0)
    dolfin_set("output destination", "terminal");
  else
    dolfin_set("output destination", "silent");

  timer1 = time(); //.restart();

  // cout << "theta: " << endl;
  // theta.vector().disp();

  W0.vector() = W.vector();
  X0.vector() = X.vector();
  rho0.vector() = rho.vector();
  S0.vector() = S.vector();
  Sdot0.vector() = Sdot.vector();
  //  P0 = P;

  mqual->meshQuality();
  cout << "FSISolver mu_min before: " << mqual->mu_min << endl;

  timer0 = time();//.restart();
  smoothMesh();

  mqual->meshQuality();
  cout << "FSISolver mu_min after: " << mqual->mu_min << endl;


  message("FSISolver timer smoother: %g", time()- timer0);
}
//-----------------------------------------------------------------------------
void NSESolver::prepareiteration()
{
  // computeX(X);
  // computeW(true);
  // computeXinc();
  // deform(Xinc);

  computeXinc();
  deform_solid(Xinc);
  computeX(X);
  computeW(true);

  P0.vector() = P.vector();

//   solutionfile << U;
//   pfile << P;

  // Compute cell mean
  ComputeStabilization(mesh(), U, nu, k,
		       delta1x, delta2x, delta1invx, *LM);
  ComputeMean(mesh(), Uc, Um, U, U0, *aM, *LM);

  computeP();
  //computeRho();
  computeS();


  reassemble = true;

  cout << "FSI prepareiteration" << endl;

  cout << "Um: " << Um.vector().norm(linf) << endl;
  cout << "Uc: " << Uc.vector().norm(linf) << endl;
  cout << "U: " << U.vector().norm(linf) << endl;
  cout << "S: " << S.vector().norm(linf) << endl;
//   U.vector().disp();
//   cout << "P: " << P.vector().norm(linf) << endl;
//   P.vector().disp();
//   cout << "delta1: " << delta1.vector().norm(linf) << endl;
//   cout << "delta2: " << delta2.vector().norm(linf) << endl;

   // std::vector<std::pair<Function*, std::string> > output;
   // std::pair<Function*, std::string> U_output(&U, "Velocity");
   // output.push_back(U_output);

   // solutionfile << output;

}
//-----------------------------------------------------------------------------
void NSESolver::postiteration()
{
}
//-----------------------------------------------------------------------------
bool NSESolver::update(real t, bool end)
{
  cout << "FSISolver::update: " << "t: " << t << " k: " << k << endl;

  // mqual->meshQuality();
  // cout << "FSISolver mu_min before: " << mqual->mu_min << endl;
  
  if(t < 30 * k)
    mu_bar = mqual->mu_min;

  if(t > 10 * k)
  {
//     if(k != 0.5*hmin/ubar)
//     {
//       cout << "increasing k to: " << k << endl;
//     }

    // Take more efficient time-step after startup
    real cfl = 1.0;
    if(mqual->mu_min < 0.75 * mu_bar)
      cfl = 0.5;
    if(mqual->mu_min < 0.5 * mu_bar)
      cfl = 0.2;
    if(mqual->mu_min < 0.4 * mu_bar)
      cfl = 0.1;

    real um = Um.vector().norm(linf);
    real wm = Wm.vector().norm(linf);
    real umax = std::max(um, wm);
    //k = std::min(cfl*hmin/umax, cfl*hmin/ubar);
    //k *= 0.2;
    //k = 2.0*hmin/ubar;
  }

//  timer0.restart();
  //smoothMesh();

  //mqual->meshQuality();
  //cout << "FSISolver mu_min after: " << mqual->mu_min << endl;
  //wfile << W;


  //message("FSISolver timer smoother: %g", time()- timer0);

  message("FSISolver timer step: %g", time()- timer1);

//   if(t > 200*k)
//     mu = 1.0e9;

  return true;
}
//-----------------------------------------------------------------------------
void NSESolver::smoothMesh()
{
  // Store mesh coordinates before smoothing
  //Xtmp.vector() = X0.vector();
  computeX(Xtmp2);
  //computeX(Xtmp);

  MeshFunction<bool> smoothed(mesh(), mesh().topology().dim());
 
  smoothed = true;
  
  //ElasticSmoother smoother(mesh());
  if(!smoother)
    smoother = new ElasticSmoother(mesh());
  if(!lsmoother)
    lsmoother = new LaplacianSmoother(mesh());


  bool did_smoothing = true;

  if(true || smooth_counter < 5)
  {
    bool reset_lsmoother = false;
    if(t == 0.0)
      reset_lsmoother = true;
    lsmoother->smooth(smoothed, solid_vertices, h0, &Wx, motionx, reset_lsmoother);
    
    Xtmp.vector() = motionx;
    Xtmp.vector() *= k;
    Xtmp.vector() += X0.vector();
    Xtmp.vector().apply();
    Xtmp.sync_ghosts();
    deform(Xtmp);
    // Wx = motionx;

    // computeX(X);
    // computeW(true);
    // computeXinc();
    // deform(Xinc);
    
    did_smoothing = true;
    
  }

  if(true || smooth_counter == 5)
  {
    int ode_max_it = dolfin_get("ODE maximum iterations");
    real ode_tol_save = dolfin_get("ODE discrete tolerance");
    dolfin_set("ODE maximum iterations", 3);
    if(mqual->mu_min < 0.5 * mu_bar)
    {
      dolfin_set("Smoother max time steps", 8);
      smoother->smooth(smoothed, solid_vertices, h0);
      did_smoothing = true;
    }
    else if(mqual->mu_min < 0.75 * mu_bar || t < 30 * k)
    {
      dolfin_set("Smoother max time steps", 2);
      smoother->smooth(smoothed, solid_vertices, h0);
      did_smoothing = true;
    }
    else
    {
      // dolfin_set("Smoother max time steps", 2);
      // smoother->smooth(smoothed, solid_vertices, h0);
      // did_smoothing = true;
    }

    if(!smoother->reset_tensor)
      smoother->J->zero();

    dolfin_set("ODE maximum iterations", ode_max_it);
    dolfin_set("ODE discrete tolerance", ode_tol_save);

    smooth_counter = 0;
  }

  smooth_counter++;

  Xtmp2.sync_ghosts();
  deform_solid(Xtmp2);

  // if(did_smoothing)
  // {
  //   computeX(X);
  //   computeW(false);
    
  //   // Revert mesh movement
  //   deform(Xtmp);

  //   //X0.vector() = Xtmp.vector();
  // }
  // else
  // {
  //   //X0.vector() = Xtmp.vector();
  // }
}
//-----------------------------------------------------------------------------
void NSESolver::solve_old()
{
}
//-----------------------------------------------------------------------------
void NSESolver::computeP()
{


  assembler->assemble(Pb, *LC, reset_tensor);
  assembler->assemble(PM, *aC, reset_tensor);
  for (uint i = 0; i < bc_con.size(); i++)
  {
    bc_con[i]->apply(PM, Pb, *aC); 
  }

  //PM.mult(Px, Presidual);
  //Presidual -= Pb;

  //Presnorm = Presidual.norm(linf);
  //cout << "Presnorm: " << Presnorm << endl;

  tic();
  //LUSolver pressure_lusolver;
  //pressure_lusolver.solve(PM, Px, Pb);
  pressure_solver.solve(PM, Px, Pb);
  cout << "pressure linear solve timer: " << toc() << endl;

  Presidual = P0.vector();
  Presidual -= P.vector();

  real relincr = 0.0;
  if(P.vector().norm(l2) > 1.0e-8)
    relincr = Presidual.norm(l2) / P.vector().norm(l2);

  incr += relincr;
  //incr += Presidual.norm(linf);
}
//-----------------------------------------------------------------------------
void NSESolver::computeRho()
{
  KrylovSolver ksolver_density(bicgstab, jacobi);

  assembler->assemble(Rhob, *LR, reset_tensor);
  assembler->assemble(RhoM, *aR, reset_tensor);
  for (uint i = 0; i < bc_rho.size(); i++)
  {
    bc_rho[i]->apply(RhoM, Rhob, *aR); 
  }

  ksolver_density.solve(RhoM, rhox, Rhob);
}
//-----------------------------------------------------------------------------
void NSESolver::computeS()
{
  tic();
  Sresidual = S.vector();

  ComputeVolInv(vol_invx);

  // cout << "vol_inv: " << endl;
  // vol_inv.vector().disp();

  Stmp.zero();
  assembler->assemble(Sdot.vector(), *LS, reset_tensor);

  S.vector() = Sdot.vector();
  S.vector() += Sdot0.vector();

  S.vector() *= 0.5*k;
  S.vector() += S0.vector();
  
  // FIXME: Zero out stress in fluid cells, why is this necessary?
  if(true)
  {
    int d = mesh().topology().dim();
    int N = mesh().numVertices();
    if(MPI::numProcesses() > 1)
      N = mesh().distdata().global_numVertices();
    int M = mesh().numCells();
    if(MPI::numProcesses() > 1)
      M = mesh().distdata().global_numCells();
    
    {
      UFC ufc(aS->form(), mesh(), aS->dofMaps());
      Cell c(mesh(), 0);
      uint local_dim = c.numEntities(0);
      uint *idx  = new uint[d * d * local_dim];
      uint *id  = new uint[d * d * local_dim];
      real *S_block = new real[d * d * local_dim];  
      
      for (CellIterator cell(mesh()); !cell.end(); ++cell)
      {
	ufc.update(*cell, mesh().distdata());
	(aS->dofMaps())[0].tabulate_dofs(idx, ufc.cell, cell->index());
	
	uint jj = 0;    
	for(int ii = 0; ii < d * d; ii++)
	{
	  S_block[jj] = 0.0;
	  id[jj++] = idx[ii];
	}
	if(!solid_cells.get(*cell))
	  S.vector().set(S_block, jj, id);
      }
      
      S.vector().apply();
      delete[] S_block;
      delete[] idx;
      delete[] id;
    }
  }
  
  message("computeS took %g seconds",toc());

  cout << "S norm: " << S.vector().norm(linf) << endl;

  Sresidual -= S.vector();

  real relincr = 0.0;
  if(S.vector().norm(l2) > 1.0e-8)
    relincr = Sresidual.norm(l2) / S.vector().norm(l2);

  incr += relincr;
  //incr += Sresidual.norm(linf);
}
//-----------------------------------------------------------------------------
void NSESolver::computeX(Function& XX)
{
  // Copy mesh coordinates into X array/function

  int d = mesh().topology().dim();
  UFC ufc(aM->form(), mesh(), aM->dofMaps());
  Cell c(mesh(), 0);
  uint local_dim = c.numEntities(0);
  uint *idx  = new uint[d * local_dim];
  uint *id  = new uint[d * local_dim];
  real *XX_block = new real[d * local_dim];  
  
  for (CellIterator cell(mesh()); !cell.end(); ++cell)
  {
    ufc.update(*cell, mesh().distdata());
    (aM->dofMaps())[0].tabulate_dofs(idx, ufc.cell, cell->index());
    
    uint ii = 0;
    uint jj = 0;    
    for(uint i = 0; i < d; i++) 
    {
      for(VertexIterator v(*cell); !v.end(); ++v, ii++) 
      {
	if (!mesh().distdata().is_ghost(v->index(), 0)) 
	{
	  XX_block[jj] = v->x()[i];
	  id[jj++] = idx[ii];
	}
	//else
	//{
	//}
      }
    }
    XX.vector().set(XX_block, jj, id);
  }
  XX.vector().apply();

  X.sync_ghosts();
  
  delete[] XX_block;
  delete[] idx;
  delete[] id;
}
//-----------------------------------------------------------------------------
void NSESolver::computeXinc()
{
  Xinc.vector() = U.vector();
  Xinc.vector() += W0.vector();
  
  Xinc.vector() *= 0.5*k;
  //Xinc.vector() *= k;
  Xinc.vector() += X0.vector();
  Xinc.vector().apply();
}
//-----------------------------------------------------------------------------
void NSESolver::computeW(bool solid)
{
  if(true || !solid)
  {
    W.vector() = X.vector();
    W.vector() -= X0.vector();
    
    cout << "W norm 0: " << W.vector().norm(linf) << endl;
    
    //W.vector() /= 0.5*k;
    //W.vector() -= W0.vector();
    W.vector() /= k;
    W.vector().apply();

    if(W.vector().norm(linf) > 4.0 * ubar)
    {
      cout << "reducing W norm: " << endl;
      W.vector() *= 4.0 * ubar / W.vector().norm(linf);
      W.vector().apply();
    }
  }

  cout << "W norm: " << W.vector().norm(linf) << endl;
  cout << "W0 norm: " << W0.vector().norm(linf) << endl;

  W.sync_ghosts();
    
}
//-----------------------------------------------------------------------------
void NSESolver::ComputeCellSize(Mesh& mesh, Vector& hvector)
{
  real *h = new real[mesh.numCells()];
  uint *rows = new uint[mesh.numCells()];

  // Compute cell size h
  hvector.init(mesh.numCells());	
  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    h[(*cell).index()] = (*cell).diameter();
    if (MPI::numProcesses() > 1)
      rows[(*cell).index()] = mesh.distdata().get_cell_global(cell->index());
    else
      rows[(*cell).index()] = (*cell).diameter();
  }

  hvector.set(h, mesh.numCells(), rows);
  hvector.apply();

  delete[] h;
  delete[] rows;
}
//-----------------------------------------------------------------------------
void NSESolver::GetMinimumCellSize(Mesh& mesh, real& hmin)
{
  // Get minimum cell diameter
  hmin = 1.0e6;
  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    if ((*cell).diameter() < hmin) hmin = (*cell).diameter();
  }

  real hmin_tmp = hmin;  
  MPI_Allreduce(&hmin_tmp, &hmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
}
//-----------------------------------------------------------------------------
void NSESolver::ComputeStabilization(Mesh& mesh, Function& w, real nu, real k, 
				     Vector& d1vector, Vector& d2vector,
				     Vector& d1invvector, Form& form)
{
  // Compute least-squares stabilizing terms: 
  //
  // if  h/nu > 1 or ny < 10^-10
  //   d1 = C1 * ( 0.5 / sqrt( 1/k^2 + |U|^2/h^2 ) )   
  //   d2 = C2 * h 
  // else 
  //   d1 = C1 * h^2  
  //   d2 = C2 * h^2  

  real C1 = 1.0;
  real C2 = 2.0;

  real kk = 0.2 * hmin / ubar;
  //real kk = k;

  UFC ufc(form.form(), mesh, form.dofMaps());
  real normw; 
  
  uint nsd = mesh.topology().dim(); 

  Cell c(mesh, 0);
  uint local_dim = c.numEntities(0);
  real *d1_block = new real[mesh.numCells()];
  real *d2_block = new real[mesh.numCells()];
  real *dinv_block = new real[mesh.numCells()];
  uint *rows2 = new uint[mesh.numCells()];
  real *w_block = new real[nsd * local_dim * mesh.numCells()];
  
  if(!indices) {
    indices = new uint[nsd * local_dim * mesh.numCells()];
    
    uint *ip = &indices[0];    
    for (CellIterator cell(mesh); !cell.end(); ++cell)
    {      
      ufc.update(*cell, mesh.distdata());            
      (form.dofMaps())[0].tabulate_dofs(ip, ufc.cell, cell->index());            
      ip += nsd * local_dim;      
    }
  }

  w.vector().get(w_block, nsd*local_dim * mesh.numCells(), indices);

  real *wp = &w_block[0];
  uint ci = 0;
  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {

    normw = 0.0;

    for(uint i =0; i < local_dim ; i++)
      normw += sqr( *(wp++) );

    normw = sqrt(normw);

    normw /= (*cell).numEntities(0);
    
    real h = (*cell).diameter();
    
    uint cid = 0;
    if ( MPI::numProcesses() > 1) 
      cid = mesh.distdata().get_cell_global((*cell).index());
    else 
      cid = (*cell).index();
    
    if (true || ((h/nu) > 1.0) || (nu < 1.0e-10) ){
      d1_block[ci] = C1 * (0.5 / sqrt( 1.0/sqr(kk) + sqr(normw/h)));
      d2_block[ci] = C2 * h;
      dinv_block[ci] = 1.0 / (d1_block[ci]);
    } else {
      d1_block[ci] = C1 * sqr(h);
      d2_block[ci] = C2 * sqr(h);
      dinv_block[ci] = 1.0 / (d1_block[ci]);
    }	
    rows2[ci++] = cid;      
  }
  
  d1vector.set(d1_block, mesh.numCells(), rows2);
  d2vector.set(d2_block, mesh.numCells(), rows2);
  d1invvector.set(dinv_block, mesh.numCells(), rows2);
  d1vector.apply();
  d2vector.apply();
  d1invvector.apply();

  delete[] d1_block;
  delete[] d2_block;
  delete[] dinv_block;
  delete[] rows2;
  delete[] w_block;
}
//-----------------------------------------------------------------------------
void NSESolver::ComputeMean(Mesh& mesh, Function& vc,
			    Function& vm, Function& v, Function& v0, 
			    Form& form, Form& form2)
{
  vc.vector() = v.vector();
  vc.vector() += v0.vector();
  vc.vector() *= 0.5;

  Cell cell_tmp(mesh, 0);
  uint nsd = mesh.topology().dim(); 
  uint local_dim = cell_tmp.numEntities(0);
  UFC ufc(form.form(), mesh, form.dofMaps()); 
  real *vc_block = new real[nsd * local_dim * mesh.numCells()];
  real *w_block = new real[nsd * local_dim * mesh.numCells()];
  real *vm_block = new real[nsd*mesh.numCells()];
  real *wm_block = new real[nsd*mesh.numCells()];
  
  if(!c_indices) {
    c_indices = new uint[nsd * mesh.numCells()];

    uint *cip = &c_indices[0];
    for(CellIterator c(mesh); !c.end(); ++c) {
      ufc.update(*c, mesh.distdata());
      (form.dofMaps())[3].tabulate_dofs(cip, ufc.cell, c->index());
      
      cip += nsd;
    }
  }
  vc.vector().get(vc_block, nsd * local_dim * mesh.numCells(), indices);
  W.vector().get(w_block, nsd * local_dim * mesh.numCells(), indices);

  uint mi = 0;
  real cellmean = 0.0;  
  real cellmean_w = 0.0;  
  uint ri = 0;
  for (CellIterator c(mesh); !c.end(); ++c)
  {

    for (uint i = 0; i < nsd; i++) {
      cellmean = 0.0;
      for (VertexIterator n(*c); !n.end(); ++n)
	cellmean += vc_block[ri];
      cellmean /= c->numEntities(0);
      cellmean_w = 0.0;
      for (VertexIterator n(*c); !n.end(); ++n)
	cellmean_w += w_block[ri];
      cellmean_w /= c->numEntities(0);
      ri++;
      vm_block[mi] = cellmean;
      wm_block[mi] = cellmean_w;
      mi++;
    }    
  }

  vm.vector().set(vm_block,nsd * mesh.numCells(), c_indices);
  vm.vector().apply();
  Wm.vector().set(wm_block,nsd * mesh.numCells(), c_indices);
  Wm.vector().apply();
  
  delete[] vc_block;
  delete[] vm_block;
  delete[] wm_block;
  delete[] w_block;

}
//-----------------------------------------------------------------------------
void NSESolver::SetInitialVelocity(Vector& xvel)
{
//   // Function for setting initial velocity, 
//   // given as a function of the coordinates (x,y,z).
//   //
//   // This function is only temporary: initial velocity 
//   // should be possible to set in the main-file. 


// //  real x,y,z;

//   real* xvelarr = xvel.vec().array();

//   for (VertexIterator vertex(mesh); !vertex.end(); ++vertex)
//   {
//     // Get coordinates of the vertex 
// //    real x = (*vertex).coord().x;
// //    real y = (*vertex).coord().y;
// //    real z = (*vertex).coord().z;
    
//     // Specify the initial velocity using (x,y,z) 
//     //
//     // Example: 
//     // 
//     // xvel((*vertex).index()) = sin(x)*cos(y)*sqrt(z);
    
//     xvelarr[(*vertex).index()] = 0.0;
//   }

//   xvel.vec().restore(xvelarr);
}
//-----------------------------------------------------------------------------                                              
void NSESolver::ComputeVolInv(Vector& vol_invx)
{
  vol_invx.init(mesh().numCells());

  real* icvarr = new real[vol_invx.local_size()];
  
  for (CellIterator cell(mesh()); !cell.end(); ++cell) 
  {
    icvarr[cell->index()] = 1.0 / (cell->volume());
  }
  vol_invx.set(icvarr);
  vol_invx.apply();
  
  delete[] icvarr;
}
//-----------------------------------------------------------------------------
void NSESolver::ComputeTimeDerivative(Mesh& mesh, Function& w, Function& w0,
				      real k, Function& dtw)
{
  dtw.vector() = w.vector();
  dtw.vector() -= w0.vector();
  dtw.vector() *= 1.0 / k;

  dtw.vector().apply();
  dtw.sync_ghosts();
}
//-----------------------------------------------------------------------------
void NSESolver::deform(Function& XX)
{
  MeshGeometry& geometry = mesh().geometry();
  
  uint d = mesh().topology().dim();
  uint N = mesh().numVertices();
  if(MPI::numProcesses() > 1)
    N = mesh().distdata().global_numVertices();
  
  UFC ufc(aM->form(), mesh(), aM->dofMaps());
  Cell c(mesh(), 0);
  uint local_dim = c.numEntities(0);
  uint *idx  = new uint[d * local_dim];
  uint *id  = new uint[d * local_dim];
  real *XX_block = new real[d * local_dim];  
  
  // Update the mesh
  for (CellIterator cell(mesh()); !cell.end(); ++cell)
  {
    ufc.update(*cell, mesh().distdata());
    (aM->dofMaps())[0].tabulate_dofs(idx, ufc.cell, cell->index());

    XX.vector().get(XX_block, d * local_dim, idx);

    uint j = 0;
    for(VertexIterator v(*cell); !v.end(); ++v)
    {
      Vertex& vertex = *v;

      //if(solid_vertices.get(vertex))
      //{
	for(unsigned int i = 0; i < d; i++)
	{
	  geometry.x(vertex.index(), i) = XX_block[i * local_dim + j];
	}
	//}
      j++;
    }
  }

  delete[] XX_block;
  delete[] idx;
  delete[] id;

  MPI_Barrier(dolfin::MPI::DOLFIN_COMM);
}
//-----------------------------------------------------------------------------
void NSESolver::deform_solid(Function& XX)
{
  MeshGeometry& geometry = mesh().geometry();
  
  uint d = mesh().topology().dim();
  uint N = mesh().numVertices();
  if(MPI::numProcesses() > 1)
    N = mesh().distdata().global_numVertices();
  
  UFC ufc(aM->form(), mesh(), aM->dofMaps());
  Cell c(mesh(), 0);
  uint local_dim = c.numEntities(0);
  uint *idx  = new uint[d * local_dim];
  uint *id  = new uint[d * local_dim];
  real *XX_block = new real[d * local_dim];  
  
  // Update the mesh
  for (CellIterator cell(mesh()); !cell.end(); ++cell)
  {
    ufc.update(*cell, mesh().distdata());
    (aM->dofMaps())[0].tabulate_dofs(idx, ufc.cell, cell->index());

    XX.vector().get(XX_block, d * local_dim, idx);

    uint j = 0;
    for(VertexIterator v(*cell); !v.end(); ++v)
    {
      Vertex& vertex = *v;

      if(solid_vertices.get(vertex))
      {
	for(unsigned int i = 0; i < d; i++)
	{
	  geometry.x(vertex.index(), i) = XX_block[i * local_dim + j];
	}
      }
      j++;
    }
  }

  delete[] XX_block;
  delete[] idx;
  delete[] id;

  MPI_Barrier(dolfin::MPI::DOLFIN_COMM);
}
//-----------------------------------------------------------------------------
