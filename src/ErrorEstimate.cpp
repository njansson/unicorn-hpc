#include <dolfin/fem/UFC.h>
#include <dolfin/mesh/RivaraRefinement.h>
#include <algorithm>
#include <map>

#include "unicorn/ErrorEstimate.h"

#ifdef HAS_MPI
#include <mpi.h>
#endif

using namespace dolfin;
using namespace unicorn;

//-----------------------------------------------------------------------------
ErrorEstimate::ErrorEstimate(Mesh& mesh, Form* Lres, Form* Lgradphi) :
  mesh(mesh), Lres_1(0), Lres_2(Lres), Lres_3(0), Lgradphi(Lgradphi),
  e_indx(mesh.numCells()), assembler(mesh),
  Rf(mesh), eif(mesh)
{
  init(mesh, Lres_1, Lres_2, Lres_3, Lgradphi, Rf, eif);

  if(MPI::numProcesses() > 1) {
    std::set<uint> eindi;
    std::map<uint, uint> mapping;
    for(CellIterator c(mesh); !c.end(); ++c)
      eindi.insert(mesh.distdata().get_cell_global(c->index()));
    
    e_indx.init_ghosted(eindi.size(), eindi, mapping);
    eindi.clear();

  }

  if(!ParameterSystem::parameters.defined("adapt_algorithm"))
    dolfin_add("adapt_algorithm", "simple");
  if(!ParameterSystem::parameters.defined("adapt_type"))
    dolfin_add("adapt_type", "cell");
}
//-----------------------------------------------------------------------------
ErrorEstimate::ErrorEstimate(Mesh& mesh, Form* Lres_1, Form* Lres_2, 
			     Form* Lres_3, Form* Lgradphi) :
  mesh(mesh), Lres_1(Lres_1), Lres_2(Lres_2), Lres_3(Lres_3), Lgradphi(Lgradphi),
  e_indx(mesh.numCells()), assembler(mesh),
  Rf(mesh), eif(mesh)
{
  init(mesh, Lres_1, Lres_2, Lres_3, Lgradphi, Rf, eif);

  if(MPI::numProcesses() > 1) {
    std::set<uint> eindi;
    std::map<uint, uint> mapping;
    for(CellIterator c(mesh); !c.end(); ++c)
      eindi.insert(mesh.distdata().get_cell_global(c->index()));
    
    e_indx.init_ghosted(eindi.size(), eindi, mapping);
    eindi.clear();

  }

  if(!ParameterSystem::parameters.defined("adapt_algorithm"))
    dolfin_add("adapt_algorithm", "simple");
  if(!ParameterSystem::parameters.defined("adapt_type"))
    dolfin_add("adapt_type", "cell");
}
//-----------------------------------------------------------------------------
ErrorEstimate::~ErrorEstimate()
{ 
  
}
//-----------------------------------------------------------------------------
void ErrorEstimate::init(Mesh& mesh, Form* Lres_1, Form* Lres_2, Form* Lres_3, 
			 Form* Lgradphi, MeshFunction<real>& Rf, MeshFunction<real>& eif) 
{

  
  if(Lres_1 != 0)
    res.init(mesh, res_1x, *Lres_1, 0);
  if(Lres_2 != 0)
    res.init(mesh, res_2x, *Lres_2, 0);
  if(Lres_3 != 0)
    res.init(mesh, res_3x, *Lres_3, 0);

  if(Lgradphi != 0)
  {
    gradphi.init(mesh, gradphix, *Lgradphi, 0);
  }
  else
  {
    gradphix.init(mesh.numCells());
  }

  Rf.init(mesh.topology().dim());
  eif.init(mesh.topology().dim());

  Rf = 0.0;
  eif = 0.0;
}
//-----------------------------------------------------------------------------
void ErrorEstimate::ComputeError(real& error)
{
  real sum = 0.0;
  real* eindarr = new real[mesh.numCells()];
  uint* rows = new uint[mesh.numCells()];

  uint ii = 0;
  for(CellIterator cell(mesh); !cell.end(); ++cell) 
    if(MPI::numProcesses() > 1)
      rows[ii++] = mesh.distdata().get_cell_global(cell->index());
    else
      rows[ii++] = cell->index();

  e_indx.get(eindarr, mesh.numCells(), rows);
  ii = 0;
  for(CellIterator cell(mesh); !cell.end(); ++cell)
  {
    real ei = eindarr[ii++];
    sum += ei;
  }

  if(MPI::processNumber() == 0)
    dolfin_set("output destination","terminal");
  real glb_sum = 0;
  MPI_Allreduce(&sum, &glb_sum, 1, MPI_DOUBLE, MPI_SUM, dolfin::MPI::DOLFIN_COMM);  
  error = sqrt(glb_sum);  
  message("ERROR = %g\n\n\n", error);
  dolfin_set("output destination","silent");

  delete[] rows;
  delete[] eindarr;
}
//-----------------------------------------------------------------------------
void ErrorEstimate::ComputeErrorIndicator(real t, real k, real T, real w)
{
  // Fixme: Right now we assume gradient is matrix
  //  int M = mesh.numCells();

  int d = mesh.topology().dim();

  // Assemble strong residual at time t
  if(Lres_1 != 0)
    assembler.assemble(res_1x, *Lres_1, false);
  if(Lres_2 != 0)
    assembler.assemble(res_2x, *Lres_2, false);
  if(Lres_3 != 0)
    assembler.assemble(res_3x, *Lres_3, false);

  // Assemble gradient of phi at time t
  if(Lgradphi != 0)
    assembler.assemble(gradphix, *Lgradphi, false);
  
  UFC ufc(Lgradphi->form(), mesh, Lgradphi->dofMaps());

  
  UFC *ufc_1 = 0;
  UFC *ufc_2 = 0;
  real* gradphiarr = new real[d * d];
  real* res_1arr = 0;
  real* res_2arr = 0;
  real* res_3arr = 0;
  if(Lres_1 != 0) { 
    res_1arr = new real[1]; 
    ufc_1 = new UFC(Lres_1->form(), mesh, Lres_1->dofMaps());
  }
  if(Lres_2 != 0) {
    res_2arr =  new real[d];  
    ufc_2 = new UFC(Lres_2->form(), mesh, Lres_2->dofMaps());
  }
  if(Lres_3 != 0)
    res_3arr = new real[d * mesh.numCells()];

  uint *eindi = new uint[mesh.numCells()];
  real* eindarr = new real[mesh.numCells()];
  uint i = 0;
  for(CellIterator c(mesh); !c.end(); ++c)
    if(MPI::numProcesses() > 1)
      eindi[i++] = mesh.distdata().get_cell_global(c->index());
    else
      eindi[i++] = c->index();
  e_indx.get(eindarr, mesh.numCells(), eindi);


  uint *gphii = new uint[d*d];
  uint *resii = new uint[d];

  uint ii = 0;
  // Compute error indicator using end-point quadrature in time
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    ufc.update(*c, mesh.distdata());
    
    //    (Lgradphi->dofMaps())[0].tabulate_dofs(&gphii[0], ufc.cell, c->index());
    (Lgradphi->dofMaps())[0].tabulate_dofs(gphii, ufc.cell, c->index());
    gradphix.get(gradphiarr, d*d, gphii);

    Cell& cell = *c;
    
    real Si = 0.0;
    int id = cell.index();
    real hmax = 0.0;
    for (EdgeIterator edge(cell); !edge.end(); ++edge)
    {
      Edge& e = *edge;
      
      real h = e.length();
      
      hmax = std::max(hmax, h);
    }

    real h = hmax;

    //    real h = MeshQuality::myDiameter(cell);
    //    real h = c->diameter();
    
    if(Lgradphi != 0)
    {
      for (int i = 0; i < d * d; i++)
      {
	Si += gradphiarr[i] * gradphiarr[i] * cell.volume();
	//	Si += gradphiarr[i * M + id] * gradphiarr[i * M + id] * cell.volume();
      }
    }
    else
      Si = 1.0;
    
    real normR = 0.0;


    if(Lres_1 != 0)
    {
      ufc_1->update(*c, mesh.distdata());
      (Lres_1->dofMaps())[0].tabulate_dofs(resii, ufc_1->cell, c->index());
      res_1x.get(res_1arr, 1, resii);

      normR += res_1arr[0] * res_1arr[0] * cell.volume();
    }
    
    if(Lres_2 != 0) 
    {
      ufc_2->update(*c, mesh.distdata());
      (Lres_2->dofMaps())[0].tabulate_dofs(resii, ufc_2->cell, c->index());
      res_2x.get(res_2arr, d, resii);
      
      for (int i = 0; i < d; i++)
	normR += res_2arr[i] * res_2arr[i] * cell.volume();
    }
    
    if(Lres_3 != 0)
      normR += res_3arr[id] * res_3arr[id] * cell.volume();
	
    real ei = 1.0 / T * k * h*h * Si * normR * w;

    eindarr[ii] += ei;

    Rf.set(*c, normR);
    eif.set(*c, eindarr[ii++]);
  }


  e_indx.set(eindarr, mesh.numCells(), eindi);
  e_indx.apply();


  delete[] gradphiarr;
  delete[] eindarr;
  delete[] eindi;
  delete[] gphii;
  delete[] resii;
  if( res_1arr) delete[] res_1arr;
  if( res_2arr) delete[] res_2arr;
  if( res_3arr) delete[] res_3arr;
  if( ufc_1)  delete ufc_1;
  if( ufc_2)  delete ufc_2;

}
//-----------------------------------------------------------------------------
void ErrorEstimate::ComputeLargestIndicators(std::vector<int>& cells,
					 real percentage)
{
  const std::string indicator_type = dolfin_get("adapt_type");
  if(indicator_type == "eind")   
    ComputeLargestIndicators_eind(cells,percentage);
  else if(indicator_type == "cell")
    ComputeLargestIndicators_cell(cells,percentage);
  else 
    dolfin::error("Unknown indicator type");
}
//-----------------------------------------------------------------------------
void ErrorEstimate::ComputeLargestIndicators_eind(std::vector<int>& cells,
						  real percentage)
{
  int N = mesh.numCells();
  real eind, sum_e, sum_e_local, max_e, max_e_local, min_e, min_e_local;
  sum_e = sum_e_local = max_e_local = 0.0;
  min_e_local = 1e6;
  
  std::vector<std::pair<int, real> > indicators(N);

  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    int id = (*cell).index();
    std::pair<int, real> p;
    p.first = id;
    uint ci = id;    
    if(MPI::numProcesses() > 1)
      ci = mesh.distdata().get_cell_global(ci);
    e_indx.get(&eind, 1, &ci);      
    p.second = eind;    
    indicators[id] = p;
    max_e_local = std::max(max_e_local, eind);
    min_e_local = std::min(min_e_local, eind);
    sum_e_local += p.second;
  }

  less_pair comp;
  std::sort(indicators.begin(), indicators.end(), comp);

  MPI_Allreduce(&sum_e_local, &sum_e, 1, MPI_DOUBLE,
		MPI_SUM, dolfin::MPI::DOLFIN_COMM);

  MPI_Allreduce(&max_e_local, &max_e, 1, MPI_DOUBLE, 
		MPI_MAX, dolfin::MPI::DOLFIN_COMM);

  MPI_Allreduce(&min_e_local, &min_e, 1, MPI_DOUBLE, 
		MPI_MIN, dolfin::MPI::DOLFIN_COMM);

  real threshold = (percentage * 0.01 * sum_e);
  real cutoff = (max_e + min_e) / 2.0;
  real acc_local, acc;
  acc_local = acc = 0.0;

  int iter = 0;
  while ( (fabs(acc - threshold) / threshold )  > 1e-2  && (iter++) < 10)
  {
    cutoff = (max_e + min_e) / 2.0;
    acc = acc_local = 0.0;
    cells.clear();

    for (int i = 0; i < N; i++) 
    {
      std::pair<int, real> p = indicators[N - 1 - i];

      cells.push_back(p.first);
      acc_local += p.second;

      if ( p.second < cutoff )
	break;     
    }

    MPI_Allreduce(&acc_local, &acc, 1, MPI_DOUBLE, 
		  MPI_SUM, dolfin::MPI::DOLFIN_COMM);
        
    ( acc > threshold ? (min_e = cutoff ) : (max_e = cutoff));    
  }
}
//-----------------------------------------------------------------------------
void ErrorEstimate::ComputeLargestIndicators_cell(std::vector<int>& cells,
					 real percentage)
{
  int N = mesh.numCells();
  int M = std::min((int)(N), 
		   (int)((real)mesh.distdata().global_numCells() * percentage * 0.01));
  
  
  if(MPI::processNumber() == 1)
    dolfin_set("output destination","terminal");
  message("Computing largest indicators");
  message("percentage: %f", percentage);
  message("N: %d", N);
  message("M: %d", M);
  dolfin_set("output destination","silent");


  std::vector<std::pair<int, real> > indicators(N);
  real eind;
  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    int id = (*cell).index();
    std::pair<int, real> p;
    p.first = id;
    uint ci = id;    
    if(MPI::numProcesses() > 1)
      ci = mesh.distdata().get_cell_global(ci);
    e_indx.get(&eind, 1, &ci);      
    p.second = eind;    
    indicators[id] = p;
  }

  less_pair comp;
  std::sort(indicators.begin(), indicators.end(), comp);


  real *local_eind = new real[M];
  for(int i = 0; i < M; i++)
  {
    std::pair<int, real> p = indicators[N - 1 - i];
    local_eind[M - 1 - i] = p.second;
  }


  /*
   *  FIXME reduce memory usage
   *  merge only half of the recived data
   */

  uint M_max, M_tot;
  MPI_Allreduce(&M, &M_max, 1, MPI_UNSIGNED, MPI_MAX, dolfin::MPI::DOLFIN_COMM);
  MPI_Allreduce(&M, &M_tot, 1, MPI_UNSIGNED, MPI_SUM, dolfin::MPI::DOLFIN_COMM);

  double *recv_eind = new double[M_max];
  double *global_eind = new double[M_tot];
  double *work = new double[M_tot];

  //  std::vector<double> global_eind;

  MPI_Status status;
  uint src,dest;
  uint rank =  MPI::processNumber();
  uint size =  MPI::numProcesses();
  uint nm = M;
  int num_recv;
  //  global_eind.insert(global_eind.begin(), local_eind, local_eind + M);
  memcpy(global_eind, local_eind, M*sizeof(real));

  for(uint i = 1; i < size; i++) {
    src =(rank - i + size) % size;
    dest = (rank + i) % size;

    MPI_Sendrecv(local_eind, M, MPI_DOUBLE, dest, 0, 
		 recv_eind, M_max, MPI_DOUBLE, src, 0, dolfin::MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_DOUBLE,&num_recv);
    //global_eind.insert(global_eind.end(), recv_eind, recv_eind + num_recv);
    merge(recv_eind, global_eind, work, num_recv, nm);
    memcpy(global_eind, work, M_tot * sizeof(real));
    nm += num_recv;
    
  }

  //  std::sort(global_eind.begin(), global_eind.end());
  cells.clear();
  int MM = (int)((real) mesh.distdata().global_numCells() * percentage * 0.01);
  int i = 0;
  for(int j = 0; j < MM; j++) {
    if( local_eind[M - 1 - i] >= global_eind[M_tot - 1 - j] ) {
      std::pair<int, real> p = indicators[N - 1 - i];
      cells.push_back(p.first);
      if( (i++) >= std::min(N, MM)) break;    
    }
  }

  dolfin_set("output destination", "terminal");
  message("%d marked cells on cpu %d", cells.size(), MPI::processNumber());
  dolfin_set("output destination", "silent");

  
  delete[] local_eind;
  delete[] recv_eind;
  delete[] global_eind;
  delete[] work;
}
//-----------------------------------------------------------------------------
void ErrorEstimate::AdaptiveRefinement(real percentage)
{

  real error = 0.0;
  ComputeError(error);

  if(MPI::processNumber() == 0)
    dolfin_set("output destination","terminal");
  message("err: %g", error);
  message("Adaptive refinement");
  message("cells before: %d", mesh.distdata().global_numCells());
  message("vertices before: %d", mesh.distdata().global_numVertices());
  dolfin_set("output destination","silent");  
  
  std::vector<int> cells;
  ComputeLargestIndicators(cells, percentage);
  
  MeshFunction<bool> cell_refinement_marker(mesh);
  cell_refinement_marker.init(mesh.topology().dim());
    
  int M = cells.size();
    
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    cell_refinement_marker.set(c->index(), false);
  }
    
  for(int i = 0; i < M; i++)
  {
    cell_refinement_marker.set(cells[i], true);
  }

  MeshFunction<real> cell_refinement_marker_r(mesh);
  cell_refinement_marker_r.init(mesh.topology().dim());
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    cell_refinement_marker_r.set(c->index(), cell_refinement_marker(*c));
  }

  File refinefile("marked.pvd");
  refinefile << cell_refinement_marker_r;

  if(MPI::processNumber() == 0)
    dolfin_set("output destination","terminal");

  const std::string refine_type = dolfin_get("adapt_algorithm");
  if(refine_type == "rivara")
    RivaraRefinement::refine(mesh, cell_refinement_marker);
  else if(refine_type == "simple")
    mesh.refine(cell_refinement_marker, true);
  else
    dolfin::error("Unknown refinement algorithm");
  dolfin_set("output destination","silent");

  if(MPI::processNumber() == 0)
    dolfin_set("output destination","terminal");
  message("cells after: %d", mesh.distdata().global_numCells());
  message("vertices after: %d", mesh.distdata().global_numVertices());
  dolfin_set("output destination","silent"); 
}
//-----------------------------------------------------------------------------
void ErrorEstimate::merge(real *a,real *b,real *res,int an,int bn)
{
  real *ap,*bp,*rp;
  ap=a;
  bp=b;
  rp=res;

  while(ap<a+an && bp<b+bn){ 
    if(*ap <= *bp){
      *rp=*ap;
      ap++;
      rp++;
    }
    else { 
      *rp=*bp;
      rp++;
      bp++;
    }
  }
  if(ap<a+an){
    do
      *rp=*ap;
    while(++rp && ++ap<a+an);
  }
  else{
    do
      *rp=*bp;
    while(++rp && ++bp<b+bn);
  }
}
//-----------------------------------------------------------------------------
