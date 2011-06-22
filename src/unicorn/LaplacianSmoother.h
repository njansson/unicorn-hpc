#ifndef __LAPLACIAN_SMOOTHER_H
#define __LAPLACIAN_SMOOTHER_H

#include <dolfin.h>
#include <cstring>
//#include <boost/numeric/ublas/matrix.hpp>
#include <unicorn/unicorn_config.h>
#include <unicorn/EquiAffineMap.h>
#include <unicorn/MeshQuality.h>
#include <unicorn/TimeDependentPDE.h>
#include <unicorn/Project.h>
#include <dolfin/fem/UFC.h>

//#include <boost/timer.hpp>

#if HAVE_SUNPERF_H
#include <sunperf.h>
#elif HAVE_SCSL_BLAS_H
#include <scsl_blas.h>
#elif HAVE_GSL_CBLAS_H
extern "C" {
#include <gsl_cblas.h>
}
#elif HAVE_CBLAS_H
extern "C" {
#include <cblas.h>
}
#endif

#define RM(row,col,nrow) ((row) + ((nrow)*(col)))

namespace dolfin { namespace unicorn
{
//-----------------------------------------------------------------------------
  class LaplacianSmoother
  {
  public:
    
    /// Constructor
    LaplacianSmoother(Mesh& mesh);
    ~LaplacianSmoother()
    {
    };

    void smooth(MeshFunction<bool>& smoothed_cells,
		MeshFunction<bool>& masked_vertices,
		MeshFunction<real>& h0,
		Vector* node_values,
		Vector& motionx);

    void maph0(Mesh& mesh, Mesh& sub, MeshFunction<int>& cell_map,
	       MeshFunction<real>& h0, MeshFunction<real>& subh0);

    static bool onBoundary(Cell& cell);
    
    static void worstElement(Mesh& mesh, int& index,
			     MeshFunction<bool>& masked_cells);
    static void elementNhood(Mesh& mesh, Cell& element,
			     MeshFunction<bool>& elements,
			     int depth);
    static void submesh(Mesh& mesh, Mesh& sub,
			MeshFunction<bool>& smoothed_cells,
			MeshFunction<int>& old2new_vertex,
			MeshFunction<int>& old2new_cell);
    
    Mesh& mesh;
    
    // Sub domain for Dirichlet boundary condition
    class DirichletBoundary : public SubDomain
    {
    public:
      bool inside(const real* x, bool on_boundary) const
      {
	if(on_boundary)
	  return true;
	else
	  return false;
      }
    };

    class MySource : public Function
    {
    public:
    MySource(Mesh& mesh) : Function(mesh)
      {
      }
      
      void eval(real* values, const real* x) const
      {
	int d = cell().dim();

	for(int i = 0; i < d; i++)
	{
	  values[i] = 0.0;
	}
      }
    };
    
    class MyBC : public Function
    {
    public:
    MyBC(Mesh& mesh) :
      Function(mesh)
      {
      }
      
      void eval(real* values, const real* x) const
      {
	int d = mesh().topology().dim();
	
	for(int i = 0; i < d; i++)
	{
	  values[i] = 0.0;
	}
      }
    };

  };
//-----------------------------------------------------------------------------
}}

#endif
