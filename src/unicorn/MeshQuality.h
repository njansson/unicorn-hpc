#ifndef __MESH_QUALITY_H
#define __MESH_QUALITY_H

#include <dolfin/common/constants.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>

namespace dolfin { namespace unicorn
{

  class MeshQuality
  {
  public:
    
    /// Constructor
    MeshQuality(Mesh& mesh);

    bool isInverted(uint& first);
    real cellQuality(Cell& cell) const;
    real meanRatio(Cell& cell) const;
    void meshQuality();
    void disp();

    // Compute cell diameter
    static real myDiameter(Cell& cell);
    static real myDiameterAvg(Cell& cell);
    static real reduce_min_scalar(real val);
    static real reduce_max_scalar(real val);
    static real reduce_avg_scalar(real val);

    Mesh& mesh;
    MeshFunction<int> orientation;

    real mu_min;
    real mu_max;
    real mu_avg;
    real mu2_min;
    real mu2_max;
    real mu2_avg;
    real h_min;
    real h_max;
    real h_avg;
    real vol_min;
    real vol_max;

    Point bbox_min;
    Point bbox_max;

  };


}}

#endif
