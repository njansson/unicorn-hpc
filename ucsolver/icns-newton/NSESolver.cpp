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
#include "unicorn/NSEResMomentum3D.h"
#include "unicorn/NSEResContinuity3D.h"
#include "unicorn/Drag3D.h"

#include "unicorn/NSEDualMomentum3D.h"
#include "unicorn/NSEDualContinuity3D.h"
#include "unicorn/NSEDualGradient3D.h"

#include <ufc.h>
#include <dolfin/fem/UFC.h>
#include <dolfin/fem/Form.h>


using namespace dolfin;
using namespace unicorn;
//-----------------------------------------------------------------------------
NSESolver::NSESolver(Mesh& mesh, Function& U, Function& U0,
		     Function& f, Function& phi, Function& beta,
		     Array<BoundaryCondition*>& bc_mom, 
		     BoundaryCondition& bc_con,
		     real T, real nu, real ubar, TimeDependent& td,
		     std::string solver_type)
  : TimeDependentPDE(mesh, bc_mom, T), U(U), U0(U0),
    f(f), phi(phi), beta(beta),
    bc_mom(bc_mom), bc_con(bc_con),
    T(T), nu(nu), ubar(ubar),
    solver_type(solver_type), errest(0), perrest(0),
    pfile("pressure.pvd"), td(td),
    pressure_solver(gmres, amg),
    ksp_pressure(0),
    startup(true), indices(0), c_indices(0)
{
  // Declare parameters
//   add("velocity file name", "velocity.pvd");
//   add("pressure file name", "pressure.pvd");
//   add("res_m file name", "res_m.pvd");
//   add("res_c file name", "res_c.pvd");
//   add("e_ind file name", "e_ind.pvd");
  
  if(!ParameterSystem::parameters.defined("Adaptive refinement percentage"))
    dolfin_add("Adaptive refinement percentage", 5.0);
  if(!ParameterSystem::parameters.defined("PDE number of samples"))
    dolfin_add("PDE number of samples", 100);

  td.sync(&t);

  // Get the number of space dimensions of the problem 
  int nsd = mesh.topology().dim();

  message("Number of space dimensions: %d",nsd);

//   real T0 = 0.0;        // start time 
//   real t  = 0.0;        // current time
//   real s = T - t;

  // Set time step (proportional to the minimum cell diameter) 
  GetMinimumCellSize(mesh, hmin);  

  // Take very conservative time-step for startup
  //k = 0.25*0.25*hmin/ubar;
  k = 0.25*hmin/ubar;
  //k = 4.0 * hmin/ubar;
  //k = 0.5 * 0.25*hmin/ubar;
  //real k = 0.05*hmin; 
  message("nu: %f",nu);
  message("ubar: %f",ubar);
  message("hmin: %f",hmin);
  message("k: %f",k);

  fnu = new Function(mesh, nu);
  //fk = new Function(mesh, k);

  fk = new TimeStepFunction(mesh);
  fk->k = &k;

  if ( nsd == 3) 
  {
    if(solver_type == "primal")
    {
      aM = new NSEMomentum3DBilinearForm(Um, *fnu, delta1, delta2, *fk);
      LM = new NSEMomentum3DLinearForm(U, U0, Uc, Um, P,
				       *fnu, delta1, delta2, f, *fk);	 
      aC = new NSEContinuity3DBilinearForm;
      LC = new NSEContinuity3DLinearForm(U, delta1inv);
    }
    else if(solver_type == "dual")
    {
      error("Not implemented yet");
    }
  }
  else 
  {
    error("Parallel Navier-Stokes solver only implemented for 3 space dimensions.");
  }

  Uc.init(mesh, Ucx, *aM, 0);
  Um.init(mesh, Umx, *LM, 4);

  U0.init(mesh, U0x, *aM, 0);
  P.init(mesh, Px, *aC, 0);
  P0.init(mesh, P0x, *aC, 0);
  delta1.init(mesh, delta1x, *LM, 7);
  delta2.init(mesh, delta2x, *LM, 7);
  delta1inv.init(mesh, delta1invx, *LM, 7);

  Presidual.init(P.vector().local_size());

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

  //TimeDependentPDE::save(U, t);
  cout << "Saving" << endl;
  
  if(t == 0.0)
  {
    solutionfile << U;
    pfile << P;
  }
  
  while(lastsample + sampleperiod < t)
  {
    lastsample = std::min(t, lastsample + sampleperiod);
    solutionfile << U;
    pfile << P;
  }
}
//-----------------------------------------------------------------------------
void NSESolver::preparestep()
{
  //  P0.vector().vec() = P.vector().vec();
  P0.vector() = P.vector();
  //  P0 = P;
}
//-----------------------------------------------------------------------------
void NSESolver::prepareiteration()
{
//   solutionfile << U;
//   pfile << P;

  // Compute cell mean
  ComputeMean(mesh(), Uc, Um, U, U0, *aM, *LM);
  ComputeStabilization(mesh(), U, nu, k,
		       delta1x, delta2x, delta1invx, *LM);


  computeP();

  cout << "Um: " << Um.vector().norm(linf) << endl;;
  cout << "Uc: " << Uc.vector().norm(linf) << endl;;
  cout << "U: " << U.vector().norm(linf) << endl;;

//   cout << "delta1: " << endl;
//   delta1x.disp();
//   cout << "delta2: " << endl;
//   delta2x.disp();

  if(startup)
    reassemble = true;
}
//-----------------------------------------------------------------------------
void NSESolver::postiteration()
{
//   tic();
//  pressure_solver.solve(PM, Px, Pb);
//   //solveLU(PM.mat(), Px.vec(), Pb.vec(), ksp_pressure);
//   cout << "exp pressure linear solve timer: " << toc() << endl;

//   real oldincr = incr;
//   incr = std::max(oldincr, Presnorm);
}
//-----------------------------------------------------------------------------
bool NSESolver::update(real t, bool end)
{
 
  //startup = false;
  if(startup)
  {
      //if(t > 0.0)
    //if(t > 1.0e-4)
    if(t > 5 * k)
    {
      // Take more efficient time-step after startup
      k = 0.5*hmin/ubar;
      korig = k;
      startup = false;

      cout << "increasing k to: " << k << endl;

      reassemble = true;
      
      cout << "stop reassembling" << endl;
    }
  }
  //return true; // assert(it's going to be ok);
}
//-----------------------------------------------------------------------------
void NSESolver::solve_old()
{
}
//-----------------------------------------------------------------------------
void NSESolver::computeP()
{
  tic();
  //if(t == 0.0)
//   if(reset_tensor)
//   {
//     assembler->assemble(PM, *aC, reset_tensor);
//     assembler->assemble(PMdummy, *aC, reset_tensor);
//     assembler->assemble(Pb, *LC, reset_tensor);
//     bc_con.apply(PM, Pb, *aC);
//     //PM.mat().lump(PMD.vec());
//     //solveLU_init(PM.mat(), Px.vec(), Pb.vec(), ksp_pressure);
//   }
  tic();
  assembler->assemble(Pb, *LC, reset_tensor);
  cout << "pressure vector timer: " << toc() << endl;

  tic();
  if(reset_tensor)
  {
    assembler->assemble(PM, *aC, reset_tensor);
    //    pressure_solver.solve(PM, Px, Pb);
    //    solveAMG_init(PM.down_cast<PETScMatrix>(), Px.down_cast<PETScVector>(),
    //		  Pb.down_cast<PETScVector>(), ksp_pressure);
    //solveLU_init(PM.mat(), Px.vec(), Pb.vec(), ksp_pressure);
  }
  cout << "pressure matrix timer: " << toc() << endl;

  cout << "Dirichlet pressure BC" << endl;
  //for (uint i = 0; i < bc_con.size(); i++)
  tic();
  bc_con.apply(PM, Pb, *aC); 
  cout << "pressure BC timer: " << toc() << endl;

  //Px.vec() = Pb.vec();
  //Px.vec().div(PMD.vec());

  PM.mult(Px, Presidual);
  Presidual -= Pb;

  Presnorm = Presidual.norm(linf);
  cout << "Presnorm: " << Presnorm << endl;

//   real oldres = res;
//   res = std::max(oldres, Presnorm);

  tic();
  pressure_solver.solve(PM, Px, Pb);
  cout << "exp pressure linear solve timer: " << toc() << endl;

//   real oldincr = incr;
//   incr = std::max(oldincr, Presnorm);


  // Experiment

//   tic();
//   KrylovSolver solver(gmres, amg);
//   solver.solve(PM, Px, Pb);
//   //solveLU(PM.mat(), Px.vec(), Pb.vec(), ksp_pressure);
//   cout << "exp pressure linear solve timer: " << toc() << endl;


  //exit(0);



//   tic();
//   KrylovSolver solver(gmres, amg);
//   solver.solve(PM, Px, Pb);
//   cout << "pressure linear solve timer: " << toc() << endl;

  //Px.zero();
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
  
  UFC ufc(form.form(), mesh, form.dofMaps());
  real normw; 
  
  Cell c(mesh, 0);
  uint local_dim = c.numEntities(0);
  real *d1_block = new real[mesh.numCells()];
  real *d2_block = new real[mesh.numCells()];
  real *dinv_block = new real[mesh.numCells()];
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
    
    if (((h/nu) > 1.0) || (nu < 1.0e-10) ){
      d1_block[ci] =  C1 * (0.5 / sqrt( 1.0/sqr(k) + sqr(normw/h)));
      d2_block[ci] = C2 * h;
      dinv_block[ci] = 1.0 / d1_block[ci];
    } else {
      d1_block[ci] = C1 * sqr(h);
      d2_block[ci] = C2 * sqr(h);
      dinv_block[ci] = 1.0 / d1_block[ci];
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
  delete[] rows2;
  delete[] w_block;
}
//-----------------------------------------------------------------------------
void NSESolver::ComputeMean(Mesh& mesh, Function& vc,
			    Function& vm, Function& v, Function& v0, 
			    Form& form, Form& form2)
{

  Cell cell_tmp(mesh, 0);
  uint nsd = mesh.topology().dim(); 
  uint local_dim = cell_tmp.numEntities(0);
  UFC ufc(form.form(), mesh, form.dofMaps()); 
  UFC ufc2(form2.form(), mesh, form2.dofMaps()); 
  real *v_block = new real[nsd * local_dim * mesh.numCells()];
  real *v0_block = new real[nsd * local_dim * mesh.numCells()];
  real *vmean_block = new real[nsd*mesh.numCells()];
   
 if(!indices) {
   indices = new uint[3 * local_dim * mesh.numCells()];
    
   uint *ip = &indices[0];    
   for (CellIterator cell(mesh); !cell.end(); ++cell)
     {      
       ufc2.update(*cell, mesh.distdata());            
       (form2.dofMaps())[0].tabulate_dofs(ip, ufc.cell, cell->index());            
       ip += 3 * local_dim;      
     }
 }

 
 if(!c_indices) {
    c_indices = new uint[nsd * mesh.numCells()];

    uint *cip = &c_indices[0];
    for(CellIterator c(mesh); !c.end(); ++c) {
      ufc.update(*c, mesh.distdata());
      (form.dofMaps())[2].tabulate_dofs(cip, ufc.cell, c->index());
      
      cip += nsd;
    }
  }

  dolfin_assert(indices);
  v.vector().get(v_block, nsd * local_dim * mesh.numCells(), indices);
  v0.vector().get(v0_block, nsd * local_dim * mesh.numCells(), indices);

  uint mi = 0;
  real cellmean = 0.0;  
  uint ri = 0;
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    for (uint i = 0; i < nsd; i++) {
      cellmean = 0.0;
      for (VertexIterator n(*c); !n.end(); ++n, ri++)
	cellmean += 0.5 * (v_block[ri] + v0_block[ri]);
      cellmean /= c->numEntities(0);
      vmean_block[mi++] = cellmean;
    }    
  }
  vm.vector().set(vmean_block,nsd * mesh.numCells(), c_indices);
  vm.vector().apply();
  
  delete[] v_block;
  delete[] v0_block;
  delete[] vmean_block;
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
void NSESolver::ComputeVolInv(Mesh& mesh, Vector& vol_inv, Form& form)
{

  UFC ufc(form.form(), mesh, form.dofMaps()); 

  real* varr = new real[mesh.numCells()];
  uint* rows = new uint[mesh.numCells()];

  std::set<uint> indices;
  std::map<uint, uint> mapping;

  // Compute cell size h                                                                                                     
  if(MPI::numProcesses() > 1)
    vol_inv.init(mesh.distdata().global_numCells());
  else
    vol_inv.init(mesh.numCells());

  uint ii = 0;
  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    varr[cell->index()] = 1.0 / cell->volume();
    ufc.update(*cell, mesh.distdata());
    (form.dofMaps())[4].tabulate_dofs(&rows[ii], 
				      ufc.cell, cell->index());
    indices.insert(rows[ii]);
    mapping[cell->index()] = rows[ii++];
  }
  vol_inv.init_ghosted(indices.size(), indices, mapping);

  vol_inv.set(varr, mesh.numCells(), rows);
  vol_inv.apply();


  indices.clear();

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
  dtw.sync_ghosts();
}
//-----------------------------------------------------------------------------
