// Copyrightx(C) 2005 Johan Hoffman.
// Licensed under the GNU GPL Version 2.
//
// Modified by Garth N. Wells 2005.
// Modified by Anders Logg 2005-2006.
// Modified by Niclas Jansson 2008-2010.
//
// First added:  2005
// Last changed: 2010-09-13

#include <cstring>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <ostream>
#include <iomanip>

#include <dolfin.h>
#include <dolfin/fem/UFC.h>

#include "unicorn/Drag3D.h"
#include "unicorn/NSESolver.h"
#include "unicorn/NSEMomentum3D.h"
#include "unicorn/NSEContinuity3D.h"
#include "unicorn/NSEResMomentum3D.h"
#include "unicorn/NSEResContinuity3D.h"
#include "unicorn/NSEDualMomentum3D.h"
#include "unicorn/NSEDualContinuity3D.h"
#include "unicorn/NSEDualGradient3D.h"
#include "unicorn/SpaceTimeFunction.h"
#include "unicorn/AdaptiveRefinement.h"

using namespace dolfin;
using namespace unicorn;

bool NSESolver::WALL_CLOCK_LIMIT = false;

//-----------------------------------------------------------------------------
NSESolver::NSESolver(Mesh& mesh, NodeNormal& node_normal,
		     Function& f, Function& beta,
		     Array<Function*>& aero_f,
                     Array<BoundaryCondition*>& bc_mom, 
                     BoundaryCondition& bc_con,
                     real T, real nu, real ubar,
		     Checkpoint& chkp, long& w_limit,
		     TimeDependent& td,
                     std::string solver_type)
  : mesh(mesh), node_normal(node_normal),
    f(f), beta(beta), aero_f(aero_f),
    bc_mom(bc_mom), bc_con(bc_con),
    T(T), nu(nu), ubar(ubar),
    chkp(chkp), w_limit(w_limit), 
    td(td), solver_type(solver_type),
    errest(0), perrest(0),
    indices(0), c_indices(0)
{
  if(!ParameterSystem::parameters.defined("adapt_percentage"))
    dolfin_add("adapt_percentage", 5.0);
  if(!ParameterSystem::parameters.defined("adapt_project"))
    dolfin_add("adapt_project", false);
  if(!ParameterSystem::parameters.defined("adapt_projected"))
    dolfin_add("adapt_projected", false);
  if(!ParameterSystem::parameters.defined("n_samples"))
    dolfin_add("n_samples", 20);
  if(!ParameterSystem::parameters.defined("krylov_method"))
    dolfin_add("krylov_method", "gmres");
  if(!ParameterSystem::parameters.defined("krylov_pc"))
    dolfin_add("krylov_pc", "amg");

}
//-----------------------------------------------------------------------------
NSESolver::~NSESolver()
{
  clear();
}
//-----------------------------------------------------------------------------
void NSESolver::clear() 
{
  if(indices)
    delete[] indices;
  indices = 0;
  
  if(c_indices)
    delete[] c_indices;
  c_indices = 0;
  
  if(errest)
    delete errest;
  errest = 0;

  if(perrest)
    delete perrest;
  perrest = 0;
}
//-----------------------------------------------------------------------------
void NSESolver::solve()
{

  if(chkp.restart()) 
  {
    if(solver_type == "primal" && chkp.id() != 0)
      return;
    else if(solver_type == "dual" && chkp.id() != 1)
      return;
  }

  if(dolfin::MPI::processNumber() == 0)
    dolfin_set("output destination","terminal");
  if(solver_type == "primal")
    message("Starting primal solver");
  else if(solver_type == "dual")
    message("Starting dual solver");
  dolfin_set("output destination","silent");

  struct sigaction sig_param;
  sig_param.sa_handler = NSESolver::sighandler;
  sigemptyset(&sig_param.sa_mask);
  sig_param.sa_flags = SA_RESTART;
  if(sigaction(SIGALRM, &sig_param, 0) < 0)
    perror("sigaction failed");

  itimerval itv;
  WALL_CLOCK_LIMIT = false;
  if ( w_limit ) 
  {
    itv.it_value.tv_sec = w_limit;
    itv.it_value.tv_usec = itv.it_interval.tv_sec = itv.it_interval.tv_usec = 0;
    if(setitimer(ITIMER_REAL, &itv, 0) < 0)
      perror("setitimer failed");
    if(dolfin::MPI::processNumber() == 0)
      dolfin_set("output destination","terminal");
    message("Wall clock limit set to %d seconds", w_limit);
    dolfin_set("output destination","silent");
  }


  // Get the number of space dimensions of the problem 
  int nsd = mesh.topology().dim();

  real T0 = 0.0;        // start time 
  real t  = 0.0;        // current time
  if(chkp.restart())
    t = chkp.restart_time();
  real s = T - t;
  td.sync(&t);

  // Set time step (proportional to the minimum cell diameter) 
  real hmin;
  GetMinimumCellSize(mesh, hmin);  

  real k = 0.15*hmin/ubar; 

 if(dolfin::MPI::processNumber() == 0)
   dolfin_set("output destination","terminal");
 message("beta: %f", (real) dolfin_get("beta"));
 message("nu: %f", nu);
 message("ubar: %f",ubar);
 message("hmin: %f",hmin);
 message("k: %f",k);
 dolfin_set("output destination","silent");
 
  Function u; // velocity
  Function up; // primal velocity
  Function upm; // Cell mean primal velocity
  Function u0; // velocity from previous time step 
  Function uc; // velocity linearized convection 
  Function um; // cell mean velocity
  Function dtu; // Time derivative of velocity
  Function dtup; // Time derivative of primal velocity
  Function p;   // pressure
  Function pp;   // primal pressure
  Function vol_inv;  // 1 / volume of element
  Function res_m;  // momentum residual
  Function delta1, delta2; // stabilization coefficients

  Function fk(mesh, k);
  Function fnu(mesh, nu);

  Function tau_1, tau_2, normal;
  
  // Initialize the bilinear and linear forms
  Form* amom = 0;
  Form* acon = 0;
  Form* Lmom = 0;
  Form* Lcon = 0;
  Form* Lres_m = 0;
  Form* LG = 0;
  Form** MF = new Form*[aero_f.size()];
  
  if ( nsd == 3 )
  {
    if(solver_type == "primal")
    {
      amom = new NSEMomentum3DBilinearForm(um,delta1,delta2,tau_1,tau_2,beta,fk,fnu);
      Lmom = new NSEMomentum3DLinearForm(um,u0,f,p,delta1,delta2,tau_1,tau_2,beta,fk,fnu);
      acon = new NSEContinuity3DBilinearForm(delta1);
      Lcon = new NSEContinuity3DLinearForm(uc, delta1);

      for(uint i = 0; i < aero_f.size(); i++)
	MF[i] = new Drag3DFunctional(*aero_f[i], dtu, u, um, p, fnu, delta1, delta2, f);

    }
    else if(solver_type == "dual")
    {
      amom = new NSEDualMomentum3DBilinearForm(um,delta1,delta2,fk,fnu,up);
      Lmom = new NSEDualMomentum3DLinearForm(um,u0,f,p,delta1,delta2,fk,fnu,up);
      acon = new NSEDualContinuity3DBilinearForm(delta1);
      Lcon = new NSEDualContinuity3DLinearForm(uc);
      Lres_m   = new NSEResMomentum3DLinearForm(pp, up, dtup, vol_inv);
      LG = new NSEDualGradient3DLinearForm(u, vol_inv);
    }
  } 
  else
  {
    error("Parallel Navier-Stokes solver only implemented for 3 space dimensions.");
  }
  
  // Create matrices and vectors 
  Matrix Amom, Acon;
  Vector bmom, bcon;

  // Initialize vectors for velocity and pressure 
  // x0vel: velocity from previous time step 
  // xcvel: linearized velocity 
  // xvel:  current velocity 
  // pvel:  current pressure 
  Vector xvel, xdtvel, x0vel, xcvel, xmvel, xpre, xppre, xpvel,xpmvel, xdtpvel;

  Vector vol_invx;         // vol_invx: needed for the computing strong residual               
  Vector res_mx;           // res_mx: needed for storing the momentum residual                 

  Vector d1vector, d2vector;
  Vector tau_1x, tau_2x, normalx;

  // Initialize vectors for the time step residuals of 
  // the momentum and continuity equations  
  uint n = mesh.numVertices() - mesh.distdata().num_ghost(0);
  
  Vector residual_mom(nsd*n);
  Vector residual_con(n);

  // Initialize algebraic solvers   
  KrylovSolver solver_con(krylov_method(dolfin_get("krylov_method")),
			  pc_type(dolfin_get("krylov_pc")));
  KrylovSolver solver_mom(krylov_method(dolfin_get("krylov_method")));

  u.init(mesh, xvel, *amom, 1);
  dtu.init(mesh, xdtvel, *amom, 1);
  u0.init(mesh, x0vel, *Lmom, 2);
  uc.init(mesh, xcvel, *Lcon, 1);
  um.init(mesh, xmvel, *amom, 2);
  p.init(mesh, xpre, *Lmom, 4);
  delta1.init(mesh, d1vector, *amom, 3);
  delta2.init(mesh, d2vector, *amom, 4);
  tau_1.init(mesh, tau_1x, *amom, 6); 
  tau_2.init(mesh, tau_2x, *amom, 6); 
  normal.init(mesh, normalx, *amom, 6); 
  
  std::vector<Function *>func;
  std::vector<Vector *> vec;
  func.push_back(&u);
  func.push_back(&u0);
  func.push_back(&uc);
  func.push_back(&um);
  func.push_back(&p);
  func.push_back(&delta1);
  func.push_back(&delta2);


  if(solver_type == "dual")
  {
    res_m.init(mesh, res_mx, *Lres_m, 0);
    res_mx.zero();
    vol_inv.init(mesh, vol_invx, *Lres_m, 4);
    pp.init(mesh, xppre, *Lres_m, 1);
    up.init(mesh, xpvel, *Lres_m, 2);
    dtup.init(mesh, xdtpvel, *Lres_m, 3);
    xdtpvel.zero();
    upm.init(mesh, xpmvel, *amom, 6);

    func.push_back(&res_m);
    func.push_back(&pp);
    func.push_back(&up);
    func.push_back(&dtup);
    func.push_back(&upm);

    // Compute the volume inverse of an element
    ComputeVolInv(mesh, vol_invx);   

  }
  else
    ComputeTangentialVectors(mesh, tau_1x, tau_2x, normalx, *amom, node_normal);
  
  // Initialize values
  x0vel.zero();
  xcvel.zero();
  xmvel.zero();
  xvel.zero();
  xpre.zero();

  u.sync_ghosts(); // velocity
  up.sync_ghosts(); // primal velocity
  upm.sync_ghosts(); // Cell mean primal velocity
  u0.sync_ghosts(); // velocity from previous time step 
  uc.sync_ghosts(); // velocity linearized convection 
  um.sync_ghosts(); // cell mean velocity
  dtu.sync_ghosts(); // Time derivative of velocity
  dtup.sync_ghosts(); // Time derivative of primal velocity
  p.sync_ghosts();   // pressure
  pp.sync_ghosts();   // primal pressure
  vol_inv.sync_ghosts();  // 1 / volume of element
  res_m.sync_ghosts();  // momentum residual
  delta1.sync_ghosts();
  delta2.sync_ghosts(); // stabilization coefficients
  res_m.sync_ghosts();
  tau_1.sync_ghosts();
  tau_2.sync_ghosts();
  normal.sync_ghosts();

  if(solver_type == "dual")                                                                    
  {
    errest = new ErrorEstimate(mesh, Lres_m, LG);                             
    vec.push_back(&errest->e_indx);    
  }

  if(chkp.restart())
  {
    chkp.load(func);
    chkp.load(vec);
  }
  
  if(solver_type == "primal" && dolfin_get("adapt_projected"))  
  {
    std::stringstream p_filename;
    p_filename << "../scratch/projected_0" << "_" << MPI::processNumber() << ".bin" << std::ends;
    File p_file(p_filename.str());
    p_file >> u.vector();
    u.vector().apply();
    

    //    dolfin_set("adapt_projected", false);

    File pp("projected.pvd");
    pp << u;

  }



  Assembler assembler(mesh);
  
  // Initialize output files 
  std::string f_fname = "aero_f.dat";

  std::vector<std::pair<Function*, std::string> > output;
  std::pair<Function*, std::string> u_output(&u, "Velocity");
  std::pair<Function*, std::string> p_output(&p, "Pressure");
  output.push_back(u_output);
  output.push_back(p_output);

  std::ostringstream output_filename;
  output_filename << solver_type << "_solution.pvd";
  
  std::stringstream p_ufilename;
  p_ufilename << "project_u" << "_" << MPI::processNumber() << ".bin" << std::ends;

  std::stringstream p_pfilename;
  p_pfilename << "project_p" << "_" << MPI::processNumber() << ".bin" << std::ends;

  File p_ufile(p_ufilename.str());
  File p_pfile(p_pfilename.str());

  File file_solution(output_filename.str(), t);
  File file_r("residual.pvd");
  File file_ei("ei.pvd");

  std::ofstream forceFile;

  if(solver_type == "primal" && MPI::processNumber() == 0)
  {
    forceFile.open(f_fname.c_str());
    forceFile.flush();
  }

  // Compute stabilization parameters
  ComputeStabilization(mesh,u0,nu,k,d1vector,d2vector, *Lmom);

  real iter_time;
  if(dolfin::MPI::processNumber() == 0)
    dolfin_set("output destination","terminal");
  message("Assembling matrix: continuity");
  dolfin_set("output destination","silent");

  // Assembling matrices 
  assembler.assemble(Acon, *acon);
  assembler.assemble(bcon, *Lcon);
  if(dolfin::MPI::processNumber() == 0)
    dolfin_set("output destination","terminal");
  message("Assembling matrix: momentum");
  dolfin_set("output destination","silent");
  assembler.assemble(Amom, *amom);
  assembler.assemble(bmom, *Lmom);

  // Initialize time-stepping parameters
  int time_step = 0;
  int sample = 0;
  int no_samples = dolfin_get("n_samples");
  if(chkp.restart())
    sample = ceil(( (t-k) * no_samples ) / T);

  // Residual, tolerance and maxmimum number of fixed point iterations
  real residual;
  real rtol = 1.0e-2;
  int iteration;
  int max_iteration = 50;  


  // Represent primal in space-time
  SpaceTimeFunction* Up = 0;
  SpaceTimeFunction* dtUp = 0;
  SpaceTimeFunction* Pp = 0;
  
  if (solver_type == "dual")
  {
    Up = new SpaceTimeFunction(mesh, up);
    
    std::vector<std::string> primal_fnames;
    
    Up->util_fileList("velocity", no_samples, primal_fnames);
    Up->util_addFiles(primal_fnames, T);
    
    Up->eval(s);
    
    
    dtUp = new SpaceTimeFunction (mesh, dtup);
    
    std::vector<std::string> dtprimal_fnames;
    
    dtUp->util_fileList("dtvelocity", no_samples, dtprimal_fnames);
    dtUp->util_addFiles(dtprimal_fnames, T);
    
    dtUp->eval(s);

    Pp = new SpaceTimeFunction (mesh, pp);
    
    std::vector<std::string> pprimal_fnames;
    
    Pp->util_fileList("pressure", no_samples, pprimal_fnames);
    Pp->util_addFiles(pprimal_fnames, T);
    
    Pp->eval(s);
  }


  // Start time-stepping
  if(dolfin::MPI::processNumber() == 0)
    dolfin_set("output destination","terminal");
  Progress prog("Time-stepping");
  dolfin_set("output destination","silent");
  while (t<T) 
  {
    
    time_step++;
    if(dolfin::MPI::processNumber() == 0)
      dolfin_set("output destination","terminal");
    message("Time step %d",time_step);
    dolfin_set("output destination","silent");
    
    s = T - t;
    if (solver_type == "dual")
    {
      Up->eval(s);
      dtUp->eval(s);
      Pp->eval(s);
    }

    // Set current velocity to velocity at previous time step 
    x0vel = xvel;

    // Initialize residual 
    residual = 2*rtol;
    iteration = 0;

    MPI::startTimer(iter_time);
    // Fix-point iteration for non-linear problem 
    while (residual > rtol && iteration < max_iteration){

      // Set linearized velocity to current velocity 
      xcvel = xvel;

      // Compute stabilization parameters
      tic();

      ComputeStabilization(mesh,u0,nu,k,d1vector,d2vector, *Lmom);

      if(dolfin::MPI::processNumber() == 0)
        dolfin_set("output destination","terminal");
      message("Compute stab took %g seconds",toc());
      dolfin_set("output destination","silent");

      // Compute time derivative of primal velocity
      ComputeTimeDerivative(mesh, u, u0, k, dtu);

      // Compute cell mean
      ComputeMean(mesh,um,uc, *amom);
      if(dolfin::MPI::processNumber() == 0)
        dolfin_set("output destination","terminal");
      message("Assemble vector: continuity");
      dolfin_set("output destination","silent");

      // Assemble continuity vector 
      assembler.assemble(Acon, *acon, false);
      assembler.assemble(bcon, *Lcon, false);
      bc_con.apply(Acon, bcon, *acon);

      if(dolfin::MPI::processNumber() == 0)
        dolfin_set("output destination","terminal");
      message("Assemble vector: momentum");
      dolfin_set("output destination","silent");

      // Assemble momentum vector 
      MPI::startTimer();
      assembler.assemble(Amom, *amom, false);
      assembler.assemble(bmom, *Lmom, false);
      for (uint i = 0; i < bc_mom.size(); i++)
        bc_mom[i]->apply(Amom, bmom, *amom);
      if(dolfin::MPI::processNumber() == 0)
        dolfin_set("output destination","terminal");
      message("Assemble momentum took %g seconds",MPI::stopTimer());
      dolfin_set("output destination","silent");

      // Compute residual for momentum equation 
      Amom.mult(xvel, residual_mom);
      residual_mom -= bmom;
      
      // Compute residual for continuity equation 
      Acon.mult(xpre,residual_con);
      residual_con -= bcon;
      if(dolfin::MPI::processNumber() == 0)
        dolfin_set("output destination","terminal");      
      message("Momentum residual  : l2 norm = %e",residual_mom.norm());
      message("continuity residual: l2 norm = %e",residual_con.norm());
      message("Total NSE residual : l2 norm = %e",sqrt(sqr(residual_mom.norm()) + sqr(residual_con.norm())));
      dolfin_set("output destination","silent");    
  
      residual = sqrt(sqr(residual_mom.norm()) + sqr(residual_con.norm()));

      if(dolfin::MPI::processNumber() == 0)
      dolfin_set("output destination","terminal");
      message("Solve linear system: continuity");
      // Solve the linear system for the continuity equation 
      tic();      
      solver_con.solve(Acon, xpre, bcon);
      dolfin_set("output destination","silent");      
      if(dolfin::MPI::processNumber() == 0)
        dolfin_set("output destination","terminal");      
      message("Linear solve took %g seconds",toc());

      message("Solve linear system: momentum");
      dolfin_set("output destination","silent");      
      // Solve the linear system for the momentum equation  
      tic();
      if(dolfin::MPI::processNumber() == 0)
        dolfin_set("output destination","terminal");            
      solver_mom.solve(Amom, xvel, bmom);
      message("Linear solve took %g seconds",toc());
      dolfin_set("output destination","silent");    
      
      iteration++;
    }
    if(dolfin::MPI::processNumber() == 0)
      dolfin_set("output destination","terminal");      
    message("Fix-point took %g seconds",MPI::stopTimer(iter_time));
    dolfin_set("output destination","silent");      
    
    if(solver_type == "dual")
    {
      real w = 1.0;
      if(s <= T / 2.0)
      	w = 0.0;

      errest->ComputeErrorIndicator(t, k, T, w);
    }

    if(residual > rtol)
      warning("NSE fixed point iteration did not converge"); 

    // Compute and output quantity of interest (aero_f)
    if(solver_type == "primal")
    {

      real force = 0.0;
      if( MPI::processNumber() == 0) 
	forceFile << std::setprecision(9) << std::setw(9) << std::setfill('0') << t << "\t";
            
      for (uint i = 0; i < aero_f.size(); i++)
      {
	force = assembler.assemble(*MF[i]);
	if( MPI::processNumber() == 0) 
	  forceFile << force  << "\t";
      }	
      if( MPI::processNumber() == 0)  {
	forceFile << "\n"; 
	forceFile.flush();     
      }
    }


    if ( (time_step == 1 || WALL_CLOCK_LIMIT) || (t > (T-T0)*(real(sample)/real(no_samples))) ){
      if(dolfin::MPI::processNumber() == 0)
	dolfin_set("output destination","terminal");      
      chkp.write(solver_type, solver_type == "dual", t+k, mesh, func, vec); // FIXME +k (^-^)
      message("Save solution to file");
      dolfin_set("output destination","silent");      
      file_solution << output;

      if(solver_type == "dual")
      {
	file_r << errest->Rf;
	file_ei << errest->eif;
      }

      // Save primal solution
      if(solver_type == "primal")
      {

	std::stringstream number;
	number << std::setfill('0') << std::setw(6) << sample;
	
	// Save primal velocity
	std::stringstream filename;
	filename << "velocity" << number.str() <<  "_" << MPI::processNumber() << ".bin" << std::ends;

	File velxmlfile(filename.str());
	velxmlfile << u.vector();

	// Save time derivative of primal velocity
	std::stringstream dtfilename;
	dtfilename << "dtvelocity" << number.str() <<  "_" << MPI::processNumber() << ".bin"  << std::ends;
	
	File dtxmlfile(dtfilename.str());
	dtxmlfile << dtu.vector();

	// Save pressure
	std::stringstream pfilename;
	pfilename << "pressure" << number.str() << "_" << MPI::processNumber() << ".bin" <<  std::ends;
	
	File pxmlfile(pfilename.str());
	pxmlfile << p.vector();

      }
      
      sample++;
    }

    // Increase time with timestep
    t = t + k;

    // Update progress
    if(dolfin::MPI::processNumber() == 0)
      dolfin_set("output destination","terminal");      
    prog = t / T;
    dolfin_set("output destination","silent");      
  }


  if(solver_type == "primal" && dolfin_get("adapt_project")) 
  {
    p_ufile << u.vector();
    p_pfile << p.vector();        
  }

  if(solver_type == "dual")
  {
    file_r << errest->Rf;
    file_ei << errest->eif;
  }

  if(dolfin::MPI::processNumber() == 0)
    dolfin_set("output destination","terminal");      
  chkp.write(solver_type, solver_type == "dual", t-k, mesh, func, vec); 
  message("save solution to file");
  file_solution << output;
  dolfin_set("output destination","silent");      
  
  if(solver_type == "primal" && MPI::processNumber() == 0) 
    forceFile.close();

  delete amom;
  delete Lmom;
  delete acon;
  delete Lcon;

  
  if(solver_type == "primal") {
    delete[] MF;
  }
  else if(solver_type == "dual") {
    delete Up;
    delete dtUp;
    delete Pp;
    delete LG; 
    delete Lres_m; 
  }

  chkp.reset();

  if(solver_type == "dual")
  {
    real percentage = dolfin_get("adapt_percentage");
    MeshFunction<bool> cell_marker;
    errest->ComputeRefinementMarkers(percentage, cell_marker);
    
    if(dolfin_get("adapt_project"))
    {

      dolfin_set("Load balancer redistribute", false);      

      Form *primal_amom = new NSEMomentum3DBilinearForm(um,delta1,delta2,tau_1,tau_2,beta,fk,fnu);
      //      Form *primal_Lmom = new NSEMomentum3DLinearForm(um,u0,f,p,delta1,delta2,tau_1,tau_2,beta,fk,fnu);

      Function p_primal, u_primal;
      Vector xp_primal, xu_primal;
      
      //      p_primal.init(mesh, xp_primal, *primal_Lmom, 4);
      u_primal.init(mesh, xu_primal, *primal_amom, 1);
      //      p_pfile >> p_primal.vector();
      //      p_primal.sync_ghosts();
      p_ufile >> u_primal.vector();
      u_primal.sync_ghosts();

      std::vector<AdaptiveRefinement::project_func>  pf;
      //      AdaptiveRefinement::form_tuple p_form(primal_Lmom, 4);
      //      AdaptiveRefinement::project_func p_project(&p_primal, p_form);

      AdaptiveRefinement::form_tuple u_form(primal_amom, 1);
      AdaptiveRefinement::project_func u_project(&u_primal, u_form);

      //      pf.push_back(p_project);    
      pf.push_back(u_project);    

      
      AdaptiveRefinement::refine_and_project(mesh, pf, cell_marker);
      dolfin_set("adapt_projected", true);

      //      delete primal_Lmom;
      //      delete primal_amom;
    }
    else
      AdaptiveRefinement::refine(mesh, cell_marker);
  }



  
}
//-----------------------------------------------------------------------------
void NSESolver::ComputeCellSize(Mesh& mesh, Vector& hvector)
{  
  //real* harr = hvector.down_cast<PETScVector>().array();
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
  hmin = 1.0e10;
  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    if ((*cell).diameter() < hmin) hmin = (*cell).diameter();
  }
  
  real hmin_tmp = hmin;
  
  MPI_Allreduce(&hmin_tmp, &hmin, 1, MPI_DOUBLE, 
		MPI_MIN, dolfin::MPI::DOLFIN_COMM);

}
//-----------------------------------------------------------------------------
void NSESolver::ComputeStabilization(Mesh& mesh, Function& w, real nu, real k, 
				     Vector& d1vector, Vector& d2vector,
				     Form& form)
{
  // Compute least-squares stabilizing terms: 
  //
  // if  h/nu > 1 or ny < 10^-10
  //   d1 = C1 * ( 0.5 / sqrt( 1/k^2 + |U|^2/h^2 ) )   
  //   d2 = C2 * h 
  // else 
  //   d1 = C1 * h^2  
  //   d2 = C2 * h^2  

  real C1 = 4.0;   
  real C2 = 2.0;   

  UFC ufc(form.form(), mesh, form.dofMaps());
  real normw; 

  Cell c(mesh, 0);
  uint local_dim = c.numEntities(0);
  real *d1_block = new real[mesh.numCells()];
  real *d2_block = new real[mesh.numCells()];
  uint *rows2 = new uint[mesh.numCells()];
  real *w_block = new real[3 * local_dim * mesh.numCells()];

  if(!indices) {
    indices = new uint[3 * local_dim * mesh.numCells()];

    uint *ip = &indices[0];
  
    for (CellIterator cell(mesh); !cell.end(); ++cell)
    {
      
      ufc.update(*cell, mesh.distdata());
            
      (form.dofMaps())[0].tabulate_dofs(ip, ufc.cell, cell->index());
            
      ip += 3 * local_dim;
      
    }
  }

  w.vector().get(w_block, 3*local_dim * mesh.numCells(), indices);

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
    
    if ( (nu < 1.0e-10) || ((h/nu) > 1.0) ){
      d1_block[ci] =  C1 * (0.5 / sqrt( 1.0/sqr(k) + sqr(normw/h)));
      d2_block[ci] = C2 * h;
    } else {
      d1_block[ci] = C1 * sqr(h);
      d2_block[ci] = C2 * sqr(h);
    }	
    rows2[ci++] = cid;      
  }
  
  d1vector.set(d1_block, mesh.numCells(), rows2);
  d2vector.set(d2_block, mesh.numCells(), rows2);
  d1vector.apply();
  d2vector.apply();

  delete[] d1_block;
  delete[] d2_block;
  delete[] rows2;
  delete[] w_block;

}
//-----------------------------------------------------------------------------
void NSESolver::ComputeMean(Mesh& mesh, Function& vmean, Function& v, 
			    Form& form)
{

  Cell cell_tmp(mesh, 0);
  uint nsd = mesh.topology().dim(); 
  uint local_dim = cell_tmp.numEntities(0);
  UFC ufc(form.form(), mesh, form.dofMaps()); 
  real *v_block = new real[nsd * local_dim * mesh.numCells()];
  real *vmean_block = new real[nsd*mesh.numCells()];
  
  if(!c_indices) {
    c_indices = new uint[nsd * mesh.numCells()];

    uint *cip = &c_indices[0];
    for(CellIterator c(mesh); !c.end(); ++c) {
      ufc.update(*c, mesh.distdata());
      (form.dofMaps())[2].tabulate_dofs(cip, ufc.cell, c->index());
      
      cip += nsd;
    }
  }

  v.vector().get(v_block, nsd * local_dim * mesh.numCells(), indices);

  uint mi = 0;
  real cellmean = 0.0;  
  uint ri = 0;
  for (CellIterator c(mesh); !c.end(); ++c)
  {

    for (uint i = 0; i < nsd; i++) {
      cellmean = 0.0;
      for (VertexIterator n(*c); !n.end(); ++n)
	cellmean += v_block[ri++];
      cellmean /= c->numEntities(0);
      vmean_block[mi++] = cellmean;
    }    
  }
  vmean.vector().set(vmean_block,nsd * mesh.numCells(), c_indices);
  vmean.vector().apply();
  
  delete[] v_block;
  delete[] vmean_block;

}
//-----------------------------------------------------------------------------               
void NSESolver::ComputeVolInv(Mesh& mesh, Vector& vol_inv)
{

  real* varr = new real[mesh.numCells()];
  uint* rows = new uint[mesh.numCells()];

  // Compute cell size h                                                                         
  vol_inv.init(mesh.numCells());

  uint ii = 0;
  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    varr[cell->index()] = 1.0 / cell->volume();
    rows[ii++] = mesh.distdata().get_cell_global(cell->index());
  }

  vol_inv.set(varr, mesh.numCells(), rows);
  vol_inv.apply();

  delete[] rows;
  delete[] varr;

}
//-----------------------------------------------------------------------------
void NSESolver::ComputeTimeDerivative(Mesh& mesh, Function& w, Function& w0,
                                      real k, Function& dtw)
{
  dtw.vector() = w.vector();
  dtw.vector() -= w0.vector();
  dtw.vector() *= 1.0 / k;

  dtw.vector().apply();
}
//-----------------------------------------------------------------------------
void NSESolver::ComputeTangentialVectors(Mesh& mesh,  Vector& tau_1, 
					 Vector& tau_2, Vector& normal,
					 Form& form, NodeNormal& node_normal)
{
  UFC ufc(form.form(), mesh, form.dofMaps());
  Cell c(mesh, 0);
  uint local_dim = c.numEntities(0);
  uint *idx  = new uint[3 * local_dim];
  uint *id  = new uint[3 * local_dim];
  real *tau_1_block = new real[3 * local_dim];  
  real *tau_2_block = new real[3 * local_dim];  
  real *normal_block = new real[3 * local_dim];

  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    
    ufc.update(*cell, mesh.distdata());
    
    (form.dofMaps())[1].tabulate_dofs(idx, ufc.cell, cell->index());
    
    uint ii = 0;
    uint jj = 0;    
    for(uint i = 0; i < 3; i++) 
    {
      for(VertexIterator v(*cell); !v.end(); ++v, ii++) 
      {
	if (!mesh.distdata().is_ghost(v->index(), 0)) 
	{
	  tau_1_block[jj] = node_normal.tau_1[i].get(*v);
	  tau_2_block[jj] = node_normal.tau_2[i].get(*v);
	  normal_block[jj] = node_normal.normal[i].get(*v);
	  id[jj++] = idx[ii];
	}
      }
    }

    tau_1.set(tau_1_block, jj, id);
    tau_2.set(tau_2_block, jj, id);
    normal.set(normal_block, jj, id);
  }

  tau_1.apply();
  tau_2.apply();
  normal.apply();
  delete[] tau_1_block;
  delete[] tau_2_block;
  delete[] normal_block;
  delete[] idx;
  delete[] id;

}
//-----------------------------------------------------------------------------
