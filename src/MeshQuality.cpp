#include <unicorn/MeshQuality.h>
#include <unicorn/EquiAffineMap.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Vertex.h>

#ifndef NO_UBLAS

using namespace dolfin::unicorn;
using namespace dolfin;


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
  real Fnorm = map.B.norm();
  
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
void MeshQuality::disp()
{
  cout << "Mesh quality:" << endl;
  cout << "mu_min: " << mu_min << endl;
  cout << "mu_max: " << mu_max << endl;
  cout << "mu_avg: " << mu_avg << endl;
  cout << "h_min: " << h_min << endl;
  cout << "h_max: " << h_max << endl;
  cout << "h_avg: " << h_avg << endl;
  cout << "bbox_min: " << bbox_min << endl;
  cout << "bbox_max: " << bbox_max << endl;
}
//-----------------------------------------------------------------------------

#endif
