#include "unicorn/ElasticSmoother.h"
#include "unicorn/MeshBC.h"
#include "unicorn/Elasticity2D.h"
#include "unicorn/ElasticityStress2D.h"
#include "unicorn/ElasticityJac2D.h"
#include "unicorn/Projection2D.h"
#include "unicorn/Elasticity3D.h"
#include "unicorn/ElasticityStress3D.h"
#include "unicorn/ElasticityJac3D.h"
#include "unicorn/Projection3D.h"

using namespace dolfin::unicorn;
using namespace dolfin;

//-----------------------------------------------------------------------------
ElasticSmoother::ElasticSmoother(Mesh& mesh) : mesh(mesh),
					       reset_tensor(true)
{
  cout << "ElasticSmoother ctor" << endl;
  if(!ParameterSystem::parameters.defined("Smoother max time steps"))
    dolfin_add("Smoother max time steps", 100);
  dolfin_set("Krylov maximum iterations", 100000);
  J = new Matrix;

}
//-----------------------------------------------------------------------------
void ElasticSmoother::maph0(Mesh& mesh, Mesh& sub,
			    MeshFunction<int>& cell_map,
			    MeshFunction<real>& h0,
			    MeshFunction<real>& subh0)
{
  subh0.init(sub, sub.topology().dim());

  for (CellIterator c(mesh); !c.end(); ++c)
  {
    Cell& cell = *c;

    //int subid = cell_map(cell);
//     cout << "subid: " << subid << endl;
//     cout << "h0: " << h0.get(cell) << endl;

    if(cell_map.get(cell) != -1)
      subh0.set(cell_map(cell), h0.get(cell));
  }
}
//-----------------------------------------------------------------------------
void ElasticSmoother::smooth(MeshFunction<bool>& smoothed_cells,
			     MeshFunction<bool>& masked_vertices,
			     MeshFunction<real>& h0)
{
  cout << "elastic smooth" << endl;
  int d = mesh.topology().dim();

  Mesh& sub = mesh;

  MeshFunction<int> vertex_map;
  MeshFunction<int> cell_map;

  MeshFunction<real>& subh0 = h0;

  MeshQuality q(sub);
  q.meshQuality();

  real E = 1.0; // Young's modulus
  real nu = 0.0; // Poisson's ratio

  real lmbdaval = E * nu / ((1 + nu) * (1 - 2 * nu));
  real muval = E / (2 * (1 + nu));

  real betaval = 1.0 * E; // Damping (sqrt(E))

  real kk = 1.0 / 20.0 * q.h_min * q.mu_min;

  cout << "k0: " << kk << endl;

  cout << "lmbda: " << lmbdaval << endl;
  cout << "mu: " << muval << endl;

  MySource f(sub);
  Density rho(sub);
  Quality qual(sub, subh0);
  MyBC bcf(sub);
  InitialValue u0(sub, subh0);

  DirichletBoundary dboundary;
  DirichletBC bc0(bcf, sub, dboundary);
  MeshBC bc1(sub, dboundary, masked_vertices);
  
  Function lmbda(sub, lmbdaval);
  Function mu(sub, muval);
  Function beta(sub, betaval);
  TimeStep_Function kf(sub);

  real T = 1.0e10;

  real ode_tol = 2.0e0;

  dolfin_set("ODE method", "dg");
  dolfin_set("ODE order", 0);
  //set("ODE implicit", false);
  dolfin_set("ODE implicit", true);
  dolfin_set("ODE nonlinear solver", "newton");
  dolfin_set("ODE linear solver", "iterative");
  dolfin_set("ODE matrix-free jacobian", false);
  dolfin_set("ODE M matrix constant", false);


  dolfin_set("ODE tolerance", ode_tol);
  dolfin_set("ODE discrete tolerance", ode_tol);
  
  dolfin_set("ODE monitor convergence", true);
  
  dolfin_set("ODE fixed time step", true);
  dolfin_set("ODE initial time step", kk);
  dolfin_set("ODE maximum time step", kk);

  dolfin_set("ODE maximum iterations", 2);
  
  dolfin_set("ODE save solution", false);
  dolfin_set("ODE solution file name", "primal.py");
  dolfin_set("ODE number of samples", 100);

  Function U;
  Function U0;
  Function B;
  Function B0;
  Function icv;

  Form* a = 0;
  Form* L = 0;
  Form* aS = 0;
  Form* LS = 0;

  Array <BoundaryCondition*> bc;
  bc.push_back(&bc0);
  bc.push_back(&bc1);

  
  if(d == 2)
  {
    a = new Elasticity2DBilinearForm(mu, lmbda, beta, qual, kf);
    L = new Elasticity2DLinearForm(U, U0, B, f, mu, lmbda, qual, kf);
    aS = new ElasticityStress2DBilinearForm();
    LS = new ElasticityStress2DLinearForm(B, U, icv);
  }
  else if(d == 3)
  {
    a = new Elasticity3DBilinearForm(mu, lmbda, beta, qual, kf);
    L = new Elasticity3DLinearForm(U, U0, B, f, mu, lmbda, qual, kf);
    aS = new ElasticityStress3DBilinearForm();
    LS = new ElasticityStress3DLinearForm(B, U, icv);
  }
  
  ElasticityPDE pde(a, L, aS, LS, sub, bc, T, U, U0, B, B0, kf, icv,
		    kk, subh0);
  pde.J = J;
  pde.reset_tensor = reset_tensor;

  //Compute solution
  //pde.timer1.restart();
  pde.solve(U, U0);
  //message("ElasticSmoother timer smooth: %g", pde.timer1.elapsed());

//   File smoothed_mesh("smoothed.xml");
//   smoothed_mesh << sub;

  delete a;
  delete L ;
  delete aS ;
  delete LS ;

  reset_tensor = false;

  if(!reset_tensor)
    J->zero();
}
//-----------------------------------------------------------------------------
bool ElasticSmoother::onBoundary(Cell& cell)
{
  int d = cell.mesh().topology().dim();

  for (FacetIterator f(cell); !f.end(); ++f)
  {
    if(f->numEntities(d) == 1)
    {
      return true;
    }
  }
  return false;
}
//-----------------------------------------------------------------------------
void ElasticSmoother::worstElement(Mesh& mesh, int& index,
				   MeshFunction<bool>& masked_cells)
{
  int d = mesh.topology().dim();
  mesh.init(d - 1, d);

  MeshQuality mqual(mesh);

  real mu_min = 1.0e12;
  index = -1;

  for (CellIterator c(mesh); !c.end(); ++c)
  {
    Cell& cell = *c;

    real qual = mqual.cellQuality(cell);
    if(qual < mu_min && !onBoundary(cell) && !masked_cells.get(cell))
      //if(qual < mu_min)
    {
      index = cell.index();
      mu_min = qual;
    }
    //mu_min = std::min(mu_min, qual);
  }
}
//-----------------------------------------------------------------------------
void ElasticSmoother::elementNhood(Mesh& mesh, Cell& element,
				   MeshFunction<bool>& elements,
				   int depth)
{
  elements.set(element, true);

  if(depth == 0)
    return;

  for (CellIterator c(element); !c.end(); ++c)
  {
    Cell& cell = *c;

    //elements.set(cell, true);
    elementNhood(mesh, cell, elements, depth - 1);
  }
}
//-----------------------------------------------------------------------------
void ElasticSmoother::submesh(Mesh& mesh, Mesh& sub,
			      MeshFunction<bool>& smoothed_cells,
			      MeshFunction<int>& old2new_vertex,
			      MeshFunction<int>& old2new_cell)
{
  
  old2new_vertex.init(mesh, 0);
  old2new_cell.init(mesh, mesh.topology().dim());

  int ncells = 0;
  int nvertices = 0;
  
  // Count cells and vertices in submesh
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    Cell& cell = *c;

    if(smoothed_cells.get(cell) == true)
    {
      ncells++;
    }
  }

  for (VertexIterator n(mesh); !n.end(); ++n)
  {
    Vertex& vertex = *n;

    bool included = false;

    for (CellIterator c(vertex); !c.end(); ++c)
    {
      Cell& cell = *c;
      
      if(smoothed_cells.get(cell) == true)
      {
	included = true;
      }
    }

    if(included)
    {
      nvertices++;
    }
  }

  
  // Get cell type
  const CellType& cell_type = mesh.type();

  unsigned int current_vertex = 0;
  unsigned int current_cell = 0;

  MeshEditor editor;
  MeshDistributedData distdata;
  editor.open(sub, cell_type.cellType(),
              mesh.topology().dim(), mesh.geometry().dim());

  // Specify number of vertices and cells
  editor.initVertices(nvertices);
  editor.initCells(ncells);

  for (VertexIterator n(mesh); !n.end(); ++n)
  {
    Vertex& vertex = *n;

    bool included = false;

    for (CellIterator c(vertex); !c.end(); ++c)
    {
      Cell& cell = *c;
      
      if(smoothed_cells.get(cell) == true)
      {
	included = true;
      }
    }

    if(included)
    {
      old2new_vertex.set(vertex.index(), current_vertex);
      editor.addVertex(current_vertex, vertex.point());
      distdata.set_map(current_vertex, mesh.distdata().get_global(vertex.index(),0), 0);
      if(mesh.distdata().is_ghost(vertex.index(), 0))
      {
	distdata.set_ghost(current_vertex, 0);
	distdata.set_ghost_owner(current_vertex, mesh.distdata().get_owner(vertex),0);
      }

      current_vertex++;
    }
    else
    {
      old2new_vertex.set(vertex.index(), -1);
    }
  }

  Array<unsigned int> cell_vertices(cell_type.numEntities(0));
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    Cell& cell = *c;

    if(smoothed_cells.get(cell) == true)
    {
      int cv_idx = 0;
      for (VertexIterator n(cell); !n.end(); ++n)
      {
	int id = old2new_vertex.get(n->index());
	if(id == -1)
	{
	  cout << "broken: " << n->index() << endl;
	}
	cell_vertices[cv_idx++] = id;
      }

      old2new_cell.set(cell.index(), current_cell);
      distdata.set_map(current_cell, mesh.distdata().get_cell_global(cell.index()), 3);
      editor.addCell(current_cell++, cell_vertices);    

    }
    else
    {
      old2new_cell.set(cell.index(), -1);
    }
  }


  editor.close();
  sub.distdata() = distdata;
  sub.distdata().invalid_numbering();
  sub.renumber();


  //  File submesh("submesh.pvd");
  //  submesh << sub;
  
  ///  File submeshfile("submesh.xml");
  //  submeshfile << sub;


}
//-----------------------------------------------------------------------------
void ElasticSmoother::ElasticityPDE::deform(Mesh& mesh, real k, Function& W)
{
  cout << "deform: " << endl;

  cout << "k: " << k << endl;

  MeshGeometry& geometry = mesh.geometry();
  
  int N = mesh.numVertices();
  
  real* Warr = new real[W.vector().local_size()];
  W.vector().get(Warr);

  // Update the mesh
  for(VertexIterator n(mesh); !n.end(); ++n)
  {
    Vertex& vertex = *n;
    
    for(unsigned int i = 0; i < mesh.topology().dim(); i++)
    {
      geometry.x(vertex.index(), i) += k * (Warr[i * N + vertex.index()]);
    }
  }

  delete[] Warr;
}
//-----------------------------------------------------------------------------
void ElasticSmoother::ElasticityPDE::deform_x(Function& X)
{
  MeshGeometry& geometry = mesh().geometry();
  
  uint d = mesh().topology().dim();
  uint N = mesh().numVertices();
  if(MPI::numProcesses() > 1)
    N = mesh().distdata().global_numVertices();
  
  UFC ufc(a().form(), mesh(), a().dofMaps());
  Cell c(mesh(), 0);
  uint local_dim = c.numEntities(0);
  uint *idx  = new uint[d * local_dim];
  uint *id  = new uint[d * local_dim];
  real *X_block = new real[d * local_dim];  
  
  // Update the mesh
  for (CellIterator cell(mesh()); !cell.end(); ++cell)
  {
    ufc.update(*cell, mesh().distdata());
    (a().dofMaps())[0].tabulate_dofs(idx, ufc.cell, cell->index());

    X.vector().get(X_block, d * local_dim, idx);

    uint j = 0;
    for(VertexIterator v(*cell); !v.end(); ++v)
    {
      Vertex& vertex = *v;

      if (true || !mesh().distdata().is_ghost(v->index(), 0)) 
      {
	for(unsigned int i = 0; i < d; i++)
	{
	  geometry.x(vertex.index(), i) = X_block[i * local_dim + j];
	}
      }
      j++;
    }
  }

  delete[] X_block;
  delete[] idx;
  delete[] id;

  MPI_Barrier(dolfin::MPI::DOLFIN_COMM);
}
//-----------------------------------------------------------------------------





