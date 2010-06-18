#include <unicorn/ElasticSmoother.h>
#include <unicorn/Elasticity2D.h>
#include <unicorn/ElasticityStress2D.h>
#include <unicorn/ElasticityJac2D.h>
#include <unicorn/Projection2D.h>
#include <unicorn/Elasticity3D.h>
#include <unicorn/ElasticityStress3D.h>
#include <unicorn/ElasticityJac3D.h>
#include <unicorn/Projection3D.h>

using namespace dolfin::unicorn;
using namespace dolfin;

//-----------------------------------------------------------------------------
ElasticSmoother::ElasticSmoother(Mesh& mesh) : mesh(mesh)
{
  if(!ParameterSystem::parameters.defined("Smoother max time steps"))
    dolfin_add("Smoother max time steps", 10);
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
			     MeshFunction<bool>& masked_cells,
			     MeshFunction<real>& h0,
			     Function& W, real k)
{
  cout << "elastic smooth" << endl;
  int d = mesh.topology().dim();

  Mesh sub, subns;

  MeshFunction<int> vertex_map;
  MeshFunction<int> cell_map;

  MeshFunction<real> subh0;
  submesh(mesh, sub, smoothed_cells, vertex_map, cell_map);

  cout << "submesh created" << endl;
  //sub = mesh;

  maph0(mesh, sub, cell_map, h0, subh0);

  //  File subh0file("subh0.xml");
  //  subh0file << subh0;

  MeshQuality q(sub);
  
  q.meshQuality();

  real E = 1.0e4; // Young's modulus
  real nu = 0.0; // Poisson's ratio

  real lmbdaval = E * nu / ((1 + nu) * (1 - 2 * nu));
  real muval = E / (2 * (1 + nu));

  //real betaval = E; // Damping (sqrt(E))
  real betaval = 1.0 * E; // Damping (sqrt(E))

  //real kk = 1.0 / 2.0 * h_min;
  //real kk = 1.0e-1;
  //  real kk = 1.0 / 2.0 * h_min / sqrt(E);
  real geom_scale = q.bbox_min.distance(q.bbox_max);
  //real kk = 0.1 * q.h_min * q.mu_min;
  //real kk = 0.2 * q.h_min * q.mu_min;
  //real kk = 0.2 * h_min;
  //real kk = 0.01 * 1.0 / 8.0 * q.mu_min;
  //real kk = 0.1;
  real kk = 0.01 * q.h_min * q.mu_min;

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
  
  Function lmbda(sub, lmbdaval);
  Function mu(sub, muval);
  Function beta(sub, betaval);
  TimeStep_Function kf(sub);
  //Function kf(sub, kk);

  //real T = 1.0e100;
  //  real T = kk;
  real T = 20.0;

  //real ode_tol = 1.0e-2;
  //real ode_tol = 1.0e-1;
  real ode_tol = 1.0e-3;

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
  dolfin_set("ODE maximum iterations", 1);
  //set("ODE fixed-point damping", 1.0e-01);
  
  dolfin_set("ODE monitor convergence", true);
  
  dolfin_set("ODE fixed time step", true);
  dolfin_set("ODE initial time step", kk);
  dolfin_set("ODE maximum time step", kk);
  
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
  Form* aJac = 0;
  Form* aS = 0;
  Form* LS = 0;
  Form* aP = 0;
  Form* LP = 0;

  Array <BoundaryCondition*> bc;
  bc.push_back(&bc0);

  
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

  //Compute solution
  pde.solve(U, U0);

  MeshGeometry& geometry = mesh.geometry();


  //  real* Warr = W.vector().vec().array();
  real* Warr = new real[W.vector().local_size()];

  int N = mesh.numVertices();

  for (VertexIterator n(mesh); !n.end(); ++n)
  {
    Vertex& vertex = *n;

    if(vertex_map(vertex) != -1)
    {
      Vertex subv(sub, vertex_map(vertex));

      for(int i = 0; i < d; i++)
      {
	real lx = subv.point()[i];
	real lx0 = geometry.x(vertex.index(), i);
	real lw = (lx - lx0) / k;
	
	Warr[i * N + vertex.index()] = lw;
	
	
	//       for(int i = 0; i < d; i++)
	// 	geometry.x(vertex.index(), i) = subv.point()[i];
	
      }
    }
  }


  W.vector().set(Warr);
  W.vector().apply();
  delete[] Warr;
  
  delete a;
  delete L ;
  //delete aJac ;
  delete aS ;
  delete LS ;
  //delete aP ;
  //delete LP ;
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
  cout << "index: " << index << endl;
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
  
  dolfin_debug("Entering create submesh");

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

  dolfin_debug("submesh created"); 
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

//   cout << "X: " << endl;
//   X.vector().vec().disp();
//   cout << "U: " << endl;
//   U.vector().vec().disp();

  MeshGeometry& geometry = mesh.geometry();
  
  int N = mesh.numVertices();
  
  //u.interpolate(vxvalues);
  
  //  real* Warr = W.vector().vec().array();

  real* Warr = new real[W.vector().local_size()];
  W.vector().get(Warr);

  // Update the mesh
  for(VertexIterator n(mesh); !n.end(); ++n)
  {
    Vertex& vertex = *n;
    
    for(unsigned int i = 0; i < mesh.topology().dim(); i++)
    {
//       cout << "before p(" << i << ", " << vertex.index() << ": " <<
//    	geometry.x(vertex.index(), i) << endl;
      geometry.x(vertex.index(), i) += k * (Warr[i * N + vertex.index()]);
      //geometry.x(vertex.index(), i) = vxvalues[i * N + vertex.index()];
//       cout << "p(" << i << ", " << vertex.index() << ": " <<
//    	geometry.x(vertex.index(), i) << endl;
    }
  }

  //  W.vector().set(Warr);
  //  W.vector().apply();
  delete[] Warr;
}
//-----------------------------------------------------------------------------