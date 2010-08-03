#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Vertex.h>
#include "unicorn/EquiAffineMap.h"
#include "unicorn/MeshQuality.h"


using namespace dolfin;
using namespace dolfin::unicorn;

MeshQuality::MeshQuality(Mesh& mesh) : mesh(mesh), orientation(mesh)
{
  // Initialize orientation
  orientation.init(mesh.topology().dim());

  for (CellIterator c(mesh); !c.end(); ++c)
  {
    Cell& cell = *c;
    orientation.set(cell.index(), (uint)cell.orientation());
  }
}

bool MeshQuality::isInverted(uint& first)
{
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    if(orientation.get(c->index()) != c->orientation())
    {
      first = c->index();
      return true;
    }
  }  
  return false;
}

real MeshQuality::cellQuality(Cell& cell) const
{
  EquiAffineMap map;
  map.update(cell);
  
  real det = fabs(map.det);

  // FIXME: Verify that it's the Frobenius norm
  // Compute Frobenius norm
  //  real Fnorm = map.B.norm();
  real Fnorm = 0.0;
  for (uint i = 0; i < 3; i ++) 
    for (uint j = 0; j < 3; j++)
      Fnorm += map.B[RM(i, j, 3)] * map.B[RM(i, j, 3)];
  Fnorm = sqrt(Fnorm);
  
  int d = cell.mesh().topology().dim();

  real mu = d * pow(det, 2.0 / d) / pow(Fnorm, 2.0);
  return mu;
}

void MeshQuality::meshQuality()
{
  mu_max = 0.0;
  mu_min = 1.0;

  h_max = 0.0;
  h_min = 1.0e12;

  real mu_sum = 0.0;
  real h_sum = 0.0;

  for (CellIterator c(mesh); !c.end(); ++c)
  {
    Cell& cell = *c;
    
    real mu = cellQuality(cell);
    real h = myDiameter(cell);
    
    mu_sum += mu;
    h_sum += h;

    mu_max = std::max(mu_max, mu);
    mu_min = std::min(mu_min, mu);
    
    h_max = std::max(h_max, h);
    h_min = std::min(h_min, h);

  }
  
  int d = mesh.topology().dim();
  for(int i = 0; i < d; i++)
  {
    bbox_min[i] = 1.0e12;
    bbox_max[i] = -1.0e12;
  }
  
  for (VertexIterator vi(mesh); !vi.end(); ++vi)
  {
    const Vertex& v = *vi;
    
    Point p = v.point();
    
    for(int i = 0; i < d; i++)
    {
      bbox_min[i] = std::min(bbox_min[i], p[i]);
      bbox_max[i] = std::max(bbox_max[i], p[i]);
    }
  }
  
  mu_avg = mu_sum / mesh.numCells();
  h_avg = h_sum / mesh.numCells();

  mu_min = reduce_min_scalar(mu_min);
  mu_max = reduce_max_scalar(mu_max);
  mu_avg = reduce_avg_scalar(mu_avg);

  h_min = reduce_min_scalar(h_min);
  h_max = reduce_max_scalar(h_max);
  h_avg = reduce_avg_scalar(h_avg);

}
//-----------------------------------------------------------------------------
real MeshQuality::myDiameter(Cell& cell)
{
  real hmax = 0.0;
  for (EdgeIterator edge(cell); !edge.end(); ++edge)
  {
    Edge& e = *edge;
    
    real h = e.length();
    
    hmax = std::max(hmax, h);
  }
  return hmax;
}
//-----------------------------------------------------------------------------
real MeshQuality::myDiameterAvg(Cell& cell)
{
  return pow(cell.volume(), 1.0 / cell.dim());
}
//-----------------------------------------------------------------------------
real MeshQuality::reduce_min_scalar(real val)
{
  real val_tmp = val;
  MPI_Allreduce(&val_tmp, &val, 1, MPI_DOUBLE, 
		MPI_MIN, dolfin::MPI::DOLFIN_COMM);

  return val;
}
//-----------------------------------------------------------------------------
real MeshQuality::reduce_max_scalar(real val)
{
  real val_tmp = val;
  MPI_Allreduce(&val_tmp, &val, 1, MPI_DOUBLE, 
		MPI_MAX, dolfin::MPI::DOLFIN_COMM);

  return val;
}
//-----------------------------------------------------------------------------
real MeshQuality::reduce_avg_scalar(real val)
{
  real val_tmp = val;
  MPI_Allreduce(&val_tmp, &val, 1, MPI_DOUBLE, 
		MPI_SUM, dolfin::MPI::DOLFIN_COMM);

  return val / dolfin::MPI::numProcesses();
}
//-----------------------------------------------------------------------------
void MeshQuality::disp()
{
  if(true || dolfin::MPI::processNumber() == 0)
  {
    cout << "Mesh quality rank " << dolfin::MPI::processNumber() << ":" << endl;
    cout << "mu_min: " << mu_min << endl;
    cout << "mu_max: " << mu_max << endl;
    cout << "mu_avg: " << mu_avg << endl;
    cout << "h_min: " << h_min << endl;
    cout << "h_max: " << h_max << endl;
    cout << "h_avg: " << h_avg << endl;
    cout << "bbox_min: " << bbox_min << endl;
    cout << "bbox_max: " << bbox_max << endl;
  }
}
//-----------------------------------------------------------------------------

