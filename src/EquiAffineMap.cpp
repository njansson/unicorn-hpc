// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU GPL Version 2.
//
// Modified by Garth N. Wells, 2006.
// Modified by Niclas Jansson, 2010.
//
// First added:  2005-05-17
// Last changed: 2010-06-13

#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Face.h>
#include "unicorn/EquiAffineMap.h"

using namespace dolfin;
using namespace dolfin::unicorn;

typedef std::pair<unsigned int, unsigned int> pr;

//-----------------------------------------------------------------------------
EquiAffineMap::EquiAffineMap() : B(0), C(0)
{
  det = 0.0;
  B = new real[3*3];
  C = new real[3*3];
}
//-----------------------------------------------------------------------------
EquiAffineMap::~EquiAffineMap()
{
  if(B)
    delete[] B;
  if(C) 
    delete[] C;
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
  real p[3] ={0.0, 0.0, 0.0};
  real P[3] ={0.0, 0.0, 0.0};

  P[0] = X; P[1] = Y; P[2] = Z;

  //  B.mult(P, p);
  for (uint i = 0; i < 3; i ++) 
    for (uint j = 0; j < 3; j++)
      p[i] += B[RM(i, j, 3)] * P[j];

  p[0] += p0.x(); p[1] += p0.y(); p[2] += p0.z();

  return Point(p[0], p[1], p[2]);
}
//-----------------------------------------------------------------------------
Point EquiAffineMap::mapinv(real X, real Y, real Z) const
{
  real p[3] ={0.0, 0.0, 0.0};
  real P[3] ={0.0, 0.0, 0.0};

  p[0] = X - p0.x(); p[1] = Y - p0.y(); p[2] = Z - p0.z();
  //  C.mult(p, P);
  for (uint i = 0; i < 3; i ++) 
    for (uint j = 0; j < 3; j++)
      p[i] += C[RM(i, j, 3)] * p[j];

  return Point(P[0], P[1], P[2]);
}
//-----------------------------------------------------------------------------
void EquiAffineMap::updateTriangle(Cell& cell)
{

  dolfin_assert(cell.dim() == 2);
  
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

  B[RM(0,0,3)] = -a * p0.x() + 2 * a * p1.x() - a * p2.x();
  B[RM(1,0,3)] = -a * p0.y() + 2 * a * p1.y() - a * p2.y();
  B[RM(2,0,3)] = 0.0;
  B[RM(0,1,3)] = p2.x() - p0.x();
  B[RM(1,1,3)] = p2.y() - p0.y();
  B[RM(2,1,3)] = 0.0;
  B[RM(0,2,3)] = 0.0;
  B[RM(1,2,3)] = 0.0;
  B[RM(2,2,3)] = 0.0;

  // Compute determinant
  det = B[RM(0,0,3)] * B[RM(1,1,3)] - B[RM(0,1,3)] * B[RM(1,0,3)];
  
  // Check determinant
  if ( fabs(det) < DOLFIN_EPS )
    error("Map from reference cell is singular.");
  
  // Compute inverse of Jacobian
  C[RM(0, 0, 3)] = B[RM(1, 1, 3)] / det;
  C[RM(0, 1, 3)] = -B[RM(0, 1, 3)] / det;
  C[RM(0, 2, 3)] = 0.0;
  C[RM(1, 0, 3)] = -B[RM(1, 0, 3)] / det;
  C[RM(1, 1, 3)] = B[RM(0, 0, 3)] / det;
  C[RM(1, 2, 3)] = 0.0;
  C[RM(2, 0, 3)] = 0.0;
  C[RM(2, 1, 3)] = 0.0;
  C[RM(2, 2, 3)] = 0.0;

  // Take absolute value of determinant
  det = fabs(det);
}
//-----------------------------------------------------------------------------
void EquiAffineMap::updateTetrahedron(Cell& cell)
{
  dolfin_assert(cell.dim() == 3);
  
   // Get coordinates
   p0 = Vertex(cell.mesh(), cell.entities(0)[0]).point();
   p1 = Vertex(cell.mesh(), cell.entities(0)[1]).point();
   p2 = Vertex(cell.mesh(), cell.entities(0)[2]).point();
   p3 = Vertex(cell.mesh(), cell.entities(0)[3]).point();
  
   // Compute Jacobian of map
  real a = 1.0 / sqrt(3.0);
  real b = 1.0 / sqrt(2.0) * a;
  real c = sqrt(3.0) / sqrt(2.0);
  
  B[RM(0, 0, 3)] = -a * p0.x() + 2 * a * p1.x() - a * p2.x();
  B[RM(0, 1, 3)] = -p0.x() + p2.x();
  B[RM(0, 2, 3)] = -b * p0.x() - b * p1.x() - b * p2.x() + c * p3.x();
  B[RM(1, 0, 3)] = -a * p0.y() + 2 * a * p1.y() - a * p2.y();
  B[RM(1, 1, 3)] = -p0.y() + p2.y();
  B[RM(1, 2, 3)] = -b * p0.y() - b * p1.y() - b * p2.y() + c * p3.y();
  B[RM(2, 0, 3)] = -a * p0.z() + 2 * a * p1.z() - a * p2.z();
  B[RM(2, 1, 3)] = -p0.z() + p2.z();
  B[RM(2, 2, 3)] = -b * p0.z() - b * p1.z() - b * p2.z() + c * p3.z();

   // Compute sub-determinants
   real d00 = B[RM(1, 1, 3)] * B[RM(2, 2, 3)] - B[RM(1, 2, 3)] * B[RM(2, 1, 3)];
   real d01 = B[RM(1, 2, 3)] * B[RM(2, 0, 3)] - B[RM(1, 0, 3)] * B[RM(2, 2, 3)];
   real d02 = B[RM(1, 0, 3)] * B[RM(2, 1, 3)] - B[RM(1, 1, 3)] * B[RM(2, 0, 3)];
  
   real d10 = B[RM(0, 2, 3)] * B[RM(2, 1, 3)] - B[RM(0, 1, 3)] * B[RM(2, 2, 3)];
   real d11 = B[RM(0, 0, 3)] * B[RM(2, 2, 3)] - B[RM(0, 2, 3)] * B[RM(2, 0, 3)];
   real d12 = B[RM(0, 1, 3)] * B[RM(2, 0, 3)] - B[RM(0, 0, 3)] * B[RM(2, 1, 3)];
  
   real d20 = B[RM(0, 1, 3)] * B[RM(1, 2, 3)] - B[RM(0, 2, 3)] * B[RM(1, 1, 3)];
   real d21 = B[RM(0, 2, 3)] * B[RM(1, 0, 3)] - B[RM(0, 0, 3)] * B[RM(1, 2, 3)];
   real d22 = B[RM(0, 0, 3)] * B[RM(1, 1, 3)] - B[RM(0, 1, 3)] * B[RM(1, 0, 3)];
  
   // Compute determinant
   det = B[RM(0, 0, 3)] * d00 + B[RM(1, 0, 3)] * d10 + B[RM(2, 0, 3)] * d20;
  
   // Check determinant
   if ( fabs(det) < DOLFIN_EPS )
     error("Map from reference cell is singular.");
  
   // Compute inverse of Jacobian
   C[RM(0, 0, 3)] = d00 / det;
   C[RM(0, 1, 3)] = d10 / det;
   C[RM(0, 2, 3)] = d20 / det;
   C[RM(1, 0, 3)] = d01 / det;
   C[RM(1, 1, 3)] = d11 / det;
   C[RM(1, 2, 3)] = d21 / det;
   C[RM(2, 0, 3)] = d02 / det;
   C[RM(2, 1, 3)] = d12 / det;
   C[RM(2, 2, 3)] = d22 / det;

   // Take absolute value of determinant
   det = fabs(det);
}
//-----------------------------------------------------------------------------


