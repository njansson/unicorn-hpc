#include <dolfin.h>
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
  int d = cell.mesh().topology().dim();

  EquiAffineMap map;
  map.update(cell);
  
  real det = fabs(map.det);

  // FIXME: Verify that it's the Frobenius norm
  // Compute Frobenius norm
  //  real Fnorm = map.B.norm();
  real Fnorm = 0.0;
  for (uint i = 0; i < d; i ++) 
    for (uint j = 0; j < d; j++)
      Fnorm += map.B[RM(i, j, 3)] * map.B[RM(i, j, 3)];
  Fnorm = sqrt(Fnorm);
  
  real mu = d * pow(det, 2.0 / d) / pow(Fnorm, 2.0);
  return mu;
}

#define isqrt3   5.77350269189625797959429519858e-01        /*  1.0/sqrt(3.0)*/
#define isqrt6   4.08248290463863052509822647505e-01        /*  1.0/sqrt(6.0)*/


inline bool m_fcn_3e(double &obj, const Point x[4])
{
  double matr[9], f;
  double g;

//   double a = 3.0;
//   double b = -1;
//   double c = 2.0 / 3.0;
  double a = 1.0/3.0;
  double b = 1;
  double c = -2.0 / 3.0;

  /* Calculate M = A*inv(W). */
  f       = x[1][0] + x[0][0];
  matr[0] = x[1][0] - x[0][0];
  matr[1] = (2.0*x[2][0] - f)*isqrt3;
  matr[2] = (3.0*x[3][0] - x[2][0] - f)*isqrt6;

  f       = x[1][1] + x[0][1];
  matr[3] = x[1][1] - x[0][1];
  matr[4] = (2.0*x[2][1] - f)*isqrt3;
  matr[5] = (3.0*x[3][1] - x[2][1] - f)*isqrt6;

  f       = x[1][2] + x[0][2];
  matr[6] = x[1][2] - x[0][2];
  matr[7] = (2.0*x[2][2] - f)*isqrt3;
  matr[8] = (3.0*x[3][2] - x[2][2] - f)*isqrt6;

  /* Calculate det(M). */
  g = matr[0]*(matr[4]*matr[8] - matr[5]*matr[7]) +
    matr[1]*(matr[5]*matr[6] - matr[3]*matr[8]) +
    matr[2]*(matr[3]*matr[7] - matr[4]*matr[6]);

  g = fabs(g);

  /* Calculate norm(M). */
  f = matr[0]*matr[0] + matr[1]*matr[1] + matr[2]*matr[2] +
      matr[3]*matr[3] + matr[4]*matr[4] + matr[5]*matr[5] +
    matr[6]*matr[6] + matr[7]*matr[7] + matr[8]*matr[8];

  /* Calculate objective function. */
  obj = a * pow(f, b) * pow(g, c);
  return true;
}


real MeshQuality::meanRatio(Cell& cell) const
{
  real mu = 0.0;

  Point x[4];
  
  int i = 0;
  for (VertexIterator vi(cell); !vi.end(); ++vi)
  {
    const Vertex& v = *vi;
    
    Point p = v.point();
    x[i] = p;
    i++;
  }

  m_fcn_3e(mu, x);

  return mu;
}

void MeshQuality::meshQuality()
{
  mu_max = 0.0;
  mu_min = 1.0;
  mu2_max = 0.0;
  mu2_min = 1.0e10;

  h_max = 0.0;
  h_min = 1.0e12;

  real mu_sum = 0.0;
  real mu2_sum = 0.0;
  real h_sum = 0.0;

  for (CellIterator c(mesh); !c.end(); ++c)
  {
    Cell& cell = *c;
    
    real mu = cellQuality(cell);
    real mu2 = meanRatio(cell);
    real h = myDiameter(cell);
    
    mu_sum += mu;
    mu2_sum += mu2;
    h_sum += h;

    mu_max = std::max(mu_max, mu);
    mu_min = std::min(mu_min, mu);
    mu2_max = std::max(mu2_max, mu2);
    mu2_min = std::min(mu2_min, mu2);
    
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
  mu2_avg = mu2_sum / mesh.numCells();
  h_avg = h_sum / mesh.numCells();

  mu_min = reduce_min_scalar(mu_min);
  mu_max = reduce_max_scalar(mu_max);
  mu_avg = reduce_avg_scalar(mu_avg);

  mu2_min = reduce_min_scalar(mu2_min);
  mu2_max = reduce_max_scalar(mu2_max);
  mu2_avg = reduce_avg_scalar(mu2_avg);

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
  if(dolfin::MPI::processNumber() == 0)
  {
    cout << "Mesh quality rank " << dolfin::MPI::processNumber() << ":" << endl;
    cout << "mu_min: " << mu_min << endl;
    cout << "mu_max: " << mu_max << endl;
    cout << "mu_avg: " << mu_avg << endl;
    cout << "mu2_min: " << mu2_min << endl;
    cout << "mu2_max: " << mu2_max << endl;
    cout << "mu2_avg: " << mu2_avg << endl;
    cout << "h_min: " << h_min << endl;
    cout << "h_max: " << h_max << endl;
    cout << "h_avg: " << h_avg << endl;
    cout << "bbox_min: " << bbox_min << endl;
    cout << "bbox_max: " << bbox_max << endl;
  }
}
//-----------------------------------------------------------------------------

