// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU GPL Version 2.
//
// Modified by Garth N. Wells, 2006.
//
// First added:  2005-05-17
// Last changed: 2006-12-06

#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Face.h>
#include <unicorn/EquiAffineMap.h>

#include <boost/numeric/ublas/matrix.hpp>

using namespace dolfin;
using namespace dolfin::unicorn;

typedef std::pair<unsigned int, unsigned int> pr;

//-----------------------------------------------------------------------------
EquiAffineMap::EquiAffineMap() : B(3, 3), C(3, 3)
{
  det = 0.0;
  
}
//-----------------------------------------------------------------------------
EquiAffineMap::~EquiAffineMap()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void EquiAffineMap::update(Cell& cell)
{
  switch ( cell.dim() )
  {
  case 2:
    updateTriangle(cell);
    break;
  case 3:
    updateTetrahedron(cell);
    break;
  default:
    error("Unknown cell type for EquiAffine map.");
  }
}
//-----------------------------------------------------------------------------
Point EquiAffineMap::operator() (real X, real Y) const
{
  return map(X, Y, 0);
}
//-----------------------------------------------------------------------------
Point EquiAffineMap::operator() (real X, real Y, real Z) const
{
  return map(X, Y, Z);
}
//-----------------------------------------------------------------------------
Point EquiAffineMap::map(real X, real Y, real Z) const
{
  uBlasVector p(3), P(3);

  P[0] = X; P[1] = Y; P[2] = Z;
  B.mult(P, p);

  p[0] += p0.x(); p[1] += p0.y(); p[2] += p0.z();

  return Point(p[0], p[1], p[2]);
}
//-----------------------------------------------------------------------------
Point EquiAffineMap::mapinv(real X, real Y, real Z) const
{
  uBlasVector p(3), P(3);

  p[0] = X - p0.x(); p[1] = Y - p0.y(); p[2] = Z - p0.z();
  C.mult(p, P);

  return Point(P[0], P[1], P[2]);
}
//-----------------------------------------------------------------------------
void EquiAffineMap::updateTriangle(Cell& cell)
{
  using namespace boost::numeric::ublas;

  assert(cell.dim() == 2);
  
  // Get coordinates
  p0 = Vertex(cell.mesh(), cell.entities(0)[0]).point();
  p1 = Vertex(cell.mesh(), cell.entities(0)[1]).point();
  p2 = Vertex(cell.mesh(), cell.entities(0)[2]).point();

  // Reset matrices
  //B = zero_matrix<real>(B.size(0), B.size(0));
  //C = zero_matrix<real>(C.size(0), C.size(0));

  // Compute Jacobian of map
//   B(0, 0) = p1.x() - p0.x(); B(1, 0) = p1.y() - p0.y();
//   B(0, 1) = p2.x() - p0.x(); B(1, 1) = p2.y() - p0.y();
  real a = 1.0 / sqrt(3.0);

  B.setitem(pr(0, 0), -a * p0.x() + 2 * a * p1.x() - a * p2.x());
  B.setitem(pr(1, 0), -a * p0.y() + 2 * a * p1.y() - a * p2.y());
  B.setitem(pr(0, 1), p2.x() - p0.x());
  B.setitem(pr(1, 1), p2.y() - p0.y());

  // Compute determinant
  det = B(0, 0) * B(1, 1) - B(0, 1) * B(1, 0);
  
  // Check determinant
  if ( fabs(det) < DOLFIN_EPS )
    error("Map from reference cell is singular.");
  
  // Compute inverse of Jacobian
  C.setitem(pr(0, 0), B(1, 1) / det); C.setitem(pr(0, 1), -B(0, 1) / det);
  C.setitem(pr(1, 0), -B(1, 0) / det); C.setitem(pr(1, 1), B(0, 0) / det);

  // Take absolute value of determinant
  det = fabs(det);
}
//-----------------------------------------------------------------------------
void EquiAffineMap::updateTetrahedron(Cell& cell)
{
  assert(cell.dim() == 3);
  
   // Get coordinates
   p0 = Vertex(cell.mesh(), cell.entities(0)[0]).point();
   p1 = Vertex(cell.mesh(), cell.entities(0)[1]).point();
   p2 = Vertex(cell.mesh(), cell.entities(0)[2]).point();
   p3 = Vertex(cell.mesh(), cell.entities(0)[3]).point();
  
   // Compute Jacobian of map
  real a = 1.0 / sqrt(3.0);
  real b = 1.0 / sqrt(2.0) * a;
  real c = sqrt(3.0) / sqrt(2.0);
  
  B.setitem(pr(0, 0), -a * p0.x() + 2 * a * p1.x() - a * p2.x());
  B.setitem(pr(0, 1), -p0.x() + p2.x());
  B.setitem(pr(0, 2), -b * p0.x() - b * p1.x() - b * p2.x() + c * p3.x());
  B.setitem(pr(1, 0), -a * p0.y() + 2 * a * p1.y() - a * p2.y());
  B.setitem(pr(1, 1), -p0.y() + p2.y());
  B.setitem(pr(1, 2), -b * p0.y() - b * p1.y() - b * p2.y() + c * p3.y());
  B.setitem(pr(2, 0), -a * p0.z() + 2 * a * p1.z() - a * p2.z());
  B.setitem(pr(2, 1), -p0.z() + p2.z());
  B.setitem(pr(2, 2), -b * p0.z() - b * p1.z() - b * p2.z() + c * p3.z());

   // Compute sub-determinants
   real d00 = B(1, 1) * B(2, 2) - B(1, 2) * B(2, 1);
   real d01 = B(1, 2) * B(2, 0) - B(1, 0) * B(2, 2);
   real d02 = B(1, 0) * B(2, 1) - B(1, 1) * B(2, 0);
  
   real d10 = B(0, 2) * B(2, 1) - B(0, 1) * B(2, 2);
   real d11 = B(0, 0) * B(2, 2) - B(0, 2) * B(2, 0);
   real d12 = B(0, 1) * B(2, 0) - B(0, 0) * B(2, 1);
  
   real d20 = B(0, 1) * B(1, 2) - B(0, 2) * B(1, 1);
   real d21 = B(0, 2) * B(1, 0) - B(0, 0) * B(1, 2);
   real d22 = B(0, 0) * B(1, 1) - B(0, 1) * B(1, 0);
  
   // Compute determinant
   det = B(0, 0) * d00 + B(1, 0) * d10 + B(2, 0) * d20;
  
   // Check determinant
   if ( fabs(det) < DOLFIN_EPS )
     error("Map from reference cell is singular.");
  
   // Compute inverse of Jacobian
   C.setitem(pr(0, 0), d00 / det);
   C.setitem(pr(0, 1), d10 / det);
   C.setitem(pr(0, 2), d20 / det);
   C.setitem(pr(1, 0), d01 / det);
   C.setitem(pr(1, 1), d11 / det);
   C.setitem(pr(1, 2), d21 / det);
   C.setitem(pr(2, 0), d02 / det);
   C.setitem(pr(2, 1), d12 / det);
   C.setitem(pr(2, 2), d22 / det);

   // Take absolute value of determinant
   det = fabs(det);
}
//-----------------------------------------------------------------------------
