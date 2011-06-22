#include "unicorn/LaplacianSmoother.h"
#include "unicorn/MeshBC.h"
#include "unicorn/Laplacian2D.h"
#include "unicorn/Laplacian3D.h"
#include "unicorn/Projection2D.h"
#include "unicorn/Projection3D.h"

using namespace dolfin::unicorn;
using namespace dolfin;

//-----------------------------------------------------------------------------
LaplacianSmoother::LaplacianSmoother(Mesh& mesh) : mesh(mesh)
{
  cout << "LaplacianSmoother ctor" << endl;
}
//-----------------------------------------------------------------------------
void LaplacianSmoother::maph0(Mesh& mesh, Mesh& sub,
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
void LaplacianSmoother::smooth(MeshFunction<bool>& smoothed_cells,
			       MeshFunction<bool>& masked_vertices,
			       MeshFunction<real>& h0,
			       Vector* node_values,
			       Vector& motionx)
{
  cout << "laplacian smooth" << endl;
  int d = mesh.topology().dim();

  Mesh& sub = mesh;

  MeshFunction<int> vertex_map;
  MeshFunction<int> cell_map;

  MeshFunction<real>& subh0 = h0;

  MeshQuality q(sub);
  q.meshQuality();

  MySource f(sub);
  MyBC bcf(sub);

  DirichletBoundary dboundary;
  DirichletBC bc0(bcf, sub, dboundary);
  MeshBC bc1(sub, dboundary, masked_vertices, node_values);
  
  Form* a = 0;
  Form* L = 0;

  if(d == 2)
  {
    a = new Laplacian2DBilinearForm;
    L = new Laplacian2DLinearForm(f);
  }
  else if(d == 3)
  {
    a = new Laplacian3DBilinearForm;
    L = new Laplacian3DLinearForm(f);
  }
  
  Matrix A;
  Vector b;

  assemble(A, *a, sub);
  assemble(b, *L, sub);

  bc0.apply(A, b, *a);
  bc1.apply(A, b, *a);

  KrylovSolver ksolver(bicgstab, sor);
  ksolver.solve(A, motionx, b);

  //Compute solution
  //pde.timer1.restart();
  /// Solve
  //message("LaplacianSmoother timer smooth: %g", pde.timer1.elapsed());

  delete a;
  delete L ;

}
//-----------------------------------------------------------------------------
bool LaplacianSmoother::onBoundary(Cell& cell)
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
