#ifndef __ELASTIC_SMOOTHER_H
#define __ELASTIC_SMOOTHER_H

#include <dolfin.h>
#include <cstring>
#include <unicorn/unicorn_config.h>
#include <unicorn/EquiAffineMap.h>
#include <unicorn/MeshQuality.h>
#include <unicorn/TimeDependentPDE.h>
#include <unicorn/Project.h>
#include <dolfin/fem/UFC.h>

#if HAVE_SUNPERF_H
#include <sunperf.h>
#elif HAVE_SCSL_CBLAS_H
#include <cmplrs/cblas.h>
#elif HAVE_GSL_CBLAS_H
extern "C" {
#include <gsl_cblas.h>
}
#elif HAVE_CBLAS_H
extern "C" {
#include <cblas.h>
}
#elif HAVE_F77_BLAS
extern "C" {
  void dgemm_(char transa, char transb, int m, int  n,  int  k,
	      double  alpha,  double *a, int lda, double *b, int
	      ldb, double beta, double *c, int ldc);
}
#endif

#define RM(row,col,nrow) ((row) + ((nrow)*(col)))

namespace dolfin { namespace unicorn
{
//-----------------------------------------------------------------------------
  class TimeStep_Function : public Function
  {
  public:
  TimeStep_Function(Mesh& mesh) : Function(mesh) {}
    void eval(real* values, const real* p) const
    {
      values[0] = *k;
    }
    
    real* k;
  };
//-----------------------------------------------------------------------------
  class ElasticSmoother
  {
  public:
    
    /// Constructor
    ElasticSmoother(Mesh& mesh);
    ~ElasticSmoother()
    {
    };

    void smooth(MeshFunction<bool>& smoothed_cells,
		MeshFunction<bool>& masked_vertices,
		MeshFunction<real>& h0);

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
    
    class Density : public Function
    {
    public:
      Density(Mesh& mesh) : Function(mesh)
      {
      }
      
      void eval(real* values, const real* x) const
      {
	values[0] = 1.0;
      }
    };
    
    class Quality : public Function
    {
    public:
    Quality(Mesh& mesh, MeshFunction<real>& h0) :
      Function(mesh), q(mesh), hq(mesh), maxqual(0.0),
	minqual(1.0e12), mqual(mesh)
      {
	real p = 2.0;

	// Compute h
	q.init(mesh.topology().dim());
	hq.init(mesh.topology().dim());
      
	//MeshQuality mqual(mesh);
	//mqual = new MeshQuality(mesh);

	for (CellIterator c(mesh); !c.end(); ++c)
	{
	  Cell& cc = *c;

	  // 	EquiAffineMap map;
	  // 	map.update(cc);
	  real h = cc.diameter();
	  real hh0 = h0.get(cc);
	  //real qual = 0.25 * (h * h) / cc.volume();
	  real qq = mqual.cellQuality(cc);
	  //qual = 1.0 / pow(qual, 2.0);
	  real qual = 1.0;

	  qual = 1.0 / pow(qq, p);

	  real hqual = 1.0;
	  
	  hqual = 1.0 + pow(fabs(h - hh0), p);
	  //qual = 1.0;
	  //  	if(qq < 0.4)
	  //  	  qual = 100.0;

	  //real totqual = qual + hqual;
	  real totqual = qual;

	  maxqual = std::max(maxqual, qual);
	  minqual = std::min(minqual, qual);

	  q.set(c->index(), totqual);
	}
      }
    
      void eval(real* values, const real* x) const
      {
	//FIXME: should we have to cast?
	Cell& c = const_cast<Cell&>(cell());

	//int d = cell().dim();

	real p = 2.0;
	//real p = 1.0;

	real val = 1.0;

	real qq = mqual.cellQuality(c);
	real qual = 1.0 / pow(qq, p);
	//real qual = q.get(cell().index());
	//real scale = 1.0 * maxqual / minqual;
	real scale = 1.0;
	val = qual / scale;
 	//val = 1.0;

	values[0] = val;
	//values[0] = 1.0;
      }

      MeshFunction<real> q;
      MeshFunction<real> hq;
      real maxqual;
      real minqual;
      MeshQuality mqual;
    };

    static real initial_B(Cell& cell, MeshFunction<real>& h0, int i)
    {
      EquiAffineMap map;
      map.update(cell);
      
      real h = h0.get(cell.index());
      //real scale = pow(cell.volume(), 1.0 / cell.dim());
      //real vscale = cell.volume() / pow(cell.diameter(), cell.dim());
      //real scale = h * vscale;
      //real scale = cell.volume() / pow(cell.diameter(), cell.dim());
      //real scale = h;
      //real scale = pow(cell.volume(), 1.0 / cell.dim()) / cell.diameter() * h;
      real scale = h / sqrt((real)cell.dim());

      real Finv[3*3];
      for (uint ii = 0; ii < 3; ii++)
	for (uint jj = 0; jj < 3; jj++)
	  Finv[RM(ii,jj,3)] = map.C[RM(ii,jj,3)] * scale;
      

      real B[3*3];

#if ((HAVE_CBLAS_H || HAVE_GSL_CBLAS_H || HAVE_SCSL_CBLAS_H))
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 3, 3, 3, 1.0, 
		  &Finv[0], 3, &Finv[0], 3, 0.0, &B[0], 3);
#elif HAVE_F77_BLAS
      dgemm_('T', 'N', 3, 3, 3, 1.0, &Finv[0], 3, &Finv[0], 3, 0.0, &B[0], 3);
      error("ElasticSmoother not supported for F77 BLAS");
#endif


      int d = cell.dim();

      int ii = i % d;
      int jj = i / d;

      return B[RM(ii, jj,3)];
    }

    class InitialValue : public Function
    {
    public:
      InitialValue(Mesh& mesh, MeshFunction<real>& h0) :
	Function(mesh), h0(h0)
      {
      }
      
      void eval(real* values, const real* x) const
      {
	real val = 0.0;

	//FIXME: should we have to cast?
	Cell& c = const_cast<Cell&>(cell());

	int d = c.dim();
	int N = 2 * d + d * d;

	EquiAffineMap map;
	map.update(c);

	real h = h0.get(c.index());
	real scale = h;
	//real scale = 1.0;

	real Finv[3*3];
	for (uint ii = 0; ii < 3; ii++)
	  for (uint jj = 0; jj < 3; jj++)
	    Finv[RM(ii,jj,3)] = map.C[RM(ii,jj,3)] * scale;
      

	real B[3*3];
	memset(&B[0], 0, 3*3*sizeof(real));

#if ((HAVE_CBLAS_H || HAVE_GSL_CBLAS_H || HAVE_SCSL_CBLAS_H))
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, 
		    3, 3, 3, 1.0, &Finv[0], 3, &Finv[0], 3, 0.0, &B[0], 3);
#elif HAVE_F77_BLAS
	dgemm_('N', 'T', 3, 3, 3, 1.0, &Finv[0], 3, &Finv[0], 3, 0.0, &B[0], 3);
#endif

	for(int i = 0; i < N; i++)
	{
	  if(i < d)
	    val = x[i];
	  else if(i >= 2 * d)
	  {
	    
	    int ii = (i - 2 * d) % d;
	    int jj = (i - 2 * d) / d;
	    
	    val = B[RM(ii, jj,3)];
	  }
	  values[i] = val;
	}
      }

      real mu;
      real lmbda;
      MeshFunction<real>& h0;
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

    Matrix* J;
    bool reset_tensor;
    
    class ElasticityPDE : public TimeDependentPDE
    {
    public:
      ElasticityPDE(Form* a, Form* L,
		    Form* aS, Form* LS,
		    Mesh& mesh,
		    Array <BoundaryCondition*>& bc, real T,
		    Function& U,
		    Function& U0,
		    Function& B,
		    Function& B0,
		    TimeStep_Function& kf,
		    Function& icv,
		    real kk,
		    MeshFunction<real>& h0) :
      TimeDependentPDE(mesh, bc, T, "ElasticSmoother"),
	aS(aS), LS(LS),
	U(U), B(B), B0(B0), kf(kf), icv(icv),
	hmin(0.0), kk(kk), h0(h0),
	qual(0), num_steps(0), max_steps(dolfin_get("Smoother max time steps")),
	min_steps(0)
      {
	kf.k = &k;
	//sampleperiod = T / 10.0;
	sampleperiod = kk;
	B.init(mesh, Bx, *aS, 0);
	B0.init(mesh, B0x, *aS, 0);
	Btmp.init(B.vector().size());
	icv.init(mesh, icell_vol, *LS, 3);
	

	//     B.init(mesh, Bx, *L, 2);
	//     B0.init(mesh, B0x, *L, 2);
	
	// X is compatible with U
	X.init(mesh, Xx, *a, 0);
	//X0x.init(X.vector().size());
	X0.init(mesh, X0x, *a, 0);

	X.sync_ghosts();
	X0.sync_ghosts();

	int d = mesh.topology().dim();
	
	vxvalues = new real[d * mesh.numVertices()];
	qual = new MeshQuality(mesh);

	init(*a, *L, 0, 0, 0);
	init(U, U0);
      }

      ~ElasticityPDE()
      {
	delete [] vxvalues;
	delete qual;
      }
      
      void computeX()
      {
	X.vector() = U.vector();
	//	X.vector().copy(U.vector(), 0, 0, X.vector().size());


	X.vector() *= k;
	X.vector() += X0x;
	X.vector().apply();
      }
      
      void ComputeInvCellVolume(Vector& icell_vol)
      {
	//icell_vol.init(mesh.numCells());
	
	real* icvarr = new real[icell_vol.local_size()];
	
	for (CellIterator cell(mesh()); !cell.end(); ++cell) 
	{
	  icvarr[cell->index()] = 1.0 / (cell->volume());
	}
	icell_vol.set(icvarr);
	icell_vol.apply();

	delete[] icvarr;
      }

      void computeB()
      {
	tic();
	
	//assembler.assemble(BM, *aS);
	//BM.mat().lump(Bm.vec());

	ComputeInvCellVolume(icell_vol);

	Btmp.zero();
	assembler->assemble(Btmp, *LS);
		
	B.vector() = Btmp;
	
	B.vector() *= k;
	B.vector()+= B0.vector();
	message("computeB took %g seconds",toc());

//  	cout << "B: " << endl;
//  	B.vector().disp();

//  	cout << "B0: " << endl;
//  	B0.vector().disp();

//  	cout << "U: " << endl;
//  	U.vector().disp();

// //  	cout << "U0: " << endl;
// //  	U0.vector().disp();

//  	cout << "icv: " << endl;
//  	icv.vector().disp();

	B.vector().apply();

      }
      
      void shift()
      {
	//cout << "shift()" << endl;
	
	X0.vector() = X.vector();
	B0.vector() = B.vector();
	// X0.vector().copy(X.vector(), 0, 0,
	//			       X0.vector().size());
	//	B0.vector().copy(B.vector(), 0, 0,
	//			       B0.vector().size());
	X0.vector().apply();
	B0.vector().apply();
      }
      
      void preparestep()
      {
	//	MPI::startTimer(timer0);
// 	int d = mesh().topology().dim();
// 	int N = mesh().numVertices();
// 	int M = mesh().numCells();

// 	real* Barr = B.vector().vec().array();

// 	for (CellIterator c(mesh()); !c.end(); ++c)
// 	{
// 	  for(int i = 0; i < d * d; i++)
// 	  {
// 	    real val = initial_B(*c, h0, i);
	    
// 	    Barr[i * M + c->index()] = val;
// 	  }
// 	}

	shift();
      }
      
      void prepareiteration()
      {
	// FIXME: BC doesn't seem to always be satisfied, why not?
// 	set("output destination", "silent");
// 	for (uint i = 0; i < bc().size(); i++)
// 	  bc()[i]->apply(M, U.vector(), a());
// 	set("output destination", "terminal");
// 	cout << "Uafter: " << endl;
// 	U.vector().vec().disp();

//  	cout << "X0 0: " << endl;
//  	X0.vector().disp();

	computeX();
	computeB();
	deform_x(X);

//  	cout << "X0 1: " << endl;
//  	X0.vector().disp();

	qual->meshQuality();

	cout << "U: " << x->norm(linf) << endl;

	reassemble = true;
      }
      
      void fu(const Vector& x, Vector& dotx, real T)
	//void fu(const uBlasVector& x, uBlasVector& dotx, real T)
      {
// 	cout << "X: " << endl;
// 	X.vector().vec().disp();
// 	cout << "U: " << endl;
// 	U.vector().vec().disp();
// 	cout << "B: " << endl;
// 	B.vector().vec().disp();
	
	assembler->assemble(dotx, L(), reset_tensor);
	//     for (uint i = 0; i < bc_mom.size(); i++)
	//       bc_mom[i]->apply(M, dotx, a());
// 	set("output destination", "silent");
// 	assembler.assemble(dotx, L());
// 	for (uint i = 0; i < bc().size(); i++)
// 	  bc()[i]->apply(M, dotx, a());
// 	set("output destination", "terminal");
	
	//cout << "dtU: " << dotx.norm(linf) << endl;
	
//  	cout << "dtU: " << endl;
//  	dotx.disp();
      }
      
      void save(Function& U, real t)
      {
	/*
	cout << "Saving at t: " << t << endl;
	//     Function U_0 = U[0];

	//     cout << "U: " << endl;
	//     U.vector().disp();

	if(t == 0.0)
	{
	  solutionfile << U;
	}
	
	if(true)
	  //while(lastsample + sampleperiod < t)
	//if(lastsample + sampleperiod < t)
	{
	  cout << "saved: " << endl;
	  lastsample = std::min(t, lastsample + sampleperiod);
	  solutionfile << U;

	  //MeshQuality qual(mesh());
      
	  qual->meshQuality();
      
	  real vol_min = 1.0e12;
	  real vol_max = 0.0;
      
	  for (CellIterator c(mesh()); !c.end(); ++c)
	  {
	    Cell& cell = *c;
	
	    vol_min = std::min(vol_min, cell.volume());
	    vol_max = std::max(vol_max, cell.volume());
	  }
      
      
	  qual->disp();
	}
	*/
	//qual->disp();
      }

      void u0(GenericVector& x)
      {
	// MPI::startTimer(timer2);
	int d = mesh().topology().dim();
	int N = mesh().numVertices();
	if(MPI::numProcesses() > 1)
	  N = mesh().distdata().global_numVertices();
	int M = mesh().numCells();
	if(MPI::numProcesses() > 1)
	  M = mesh().distdata().global_numCells();

	{
	  UFC ufc(aS->form(), mesh(), aS->dofMaps());
	  Cell c(mesh(), 0);
	  uint local_dim = c.numEntities(0);
	  uint *idx  = new uint[d * d * local_dim];
	  uint *id  = new uint[d * d * local_dim];
	  real *B_block = new real[d * d * local_dim];  

	  for (CellIterator cell(mesh()); !cell.end(); ++cell)
	  {
	    ufc.update(*cell, mesh().distdata());
	    (aS->dofMaps())[0].tabulate_dofs(idx, ufc.cell, cell->index());

	    uint jj = 0;    
	    for(int ii = 0; ii < d * d; ii++)
	    {
	      real val = initial_B(*cell, h0, ii);
	      B_block[jj] = val;
	      id[jj++] = idx[ii];
	    }
	    B.vector().set(B_block, jj, id);
	  }
	  
	  B.vector().apply();
	  delete[] idx;
	  delete[] id;
	  delete[] B_block;
	}

	{
	  UFC ufc(a().form(), mesh(), a().dofMaps());
	  Cell c(mesh(), 0);
	  uint local_dim = c.numEntities(0);
	  uint *idx  = new uint[d * local_dim];
	  uint *id  = new uint[d * local_dim];
	  real *X_block = new real[d * local_dim];  
	  
	  for (CellIterator cell(mesh()); !cell.end(); ++cell)
	  {
	    ufc.update(*cell, mesh().distdata());
	    (a().dofMaps())[0].tabulate_dofs(idx, ufc.cell, cell->index());
	    
	    uint ii = 0;
	    uint jj = 0;    
	    for(uint i = 0; i < d; i++) 
	    {
	      for(VertexIterator v(*cell); !v.end(); ++v, ii++) 
	      {
		if (!mesh().distdata().is_ghost(v->index(), 0)) 
		{
		  X_block[jj] = v->x()[i];
		  id[jj++] = idx[ii];
		}
		else
		{
		}
	      }
	    }
	    X.vector().set(X_block, jj, id);
	  }
	  X.vector().apply();
	  
	  delete[] X_block;
	  delete[] idx;
	  delete[] id;
	}

	U.vector().zero();

	//message("ElasticSmoother timer u0: %g", MPI::stopTimer(timer2));
      }

      void revert()
      {
	B.vector() = B0.vector();
	B.vector().apply();
      }

      //bool update(const uBlasVector& u, real t, bool end)
      bool update(real t, bool end)
      {
	cout << "ElasticSmoother::update:" << endl;

	max_steps = dolfin_get("Smoother max time steps");

// 	cout << "U:" << endl;
// 	U.vector().disp();
// 	cout << "B:" << endl;
// 	B.vector().disp();

	//deform(mesh(), k, U);

	qual->meshQuality();
	//qual.disp();
	cout << "smoother_mu_min: " << qual->mu_min << endl;
	uint inverted_cell;
	
	if(qual->isInverted(inverted_cell)) {
	  cout << "inverted cell: " << inverted_cell << endl;
	}

	int d = mesh().topology().dim();
	int N = mesh().numVertices();

	real sstate_tol = 1.0e-9;

	real norm = dotx->norm(linf);
	cout << "dotx norm: " << norm << endl;

	real xnorm = U.vector().norm(linf);
	cout << "x norm: " << xnorm << endl;

	real hhmin = 1.0e9;

	for(CellIterator ci(mesh()); !ci.end(); ++ci)
	{
	  real h = pow(ci->volume(), 1.0 / ci->dim());
	  real qi = qual->cellQuality(*ci);

	  hhmin = std::min(hhmin, h);
	}
	hhmin = MeshQuality::reduce_min_scalar(hhmin);


	cout << "hhmin: " << hhmin << endl;

        k = 1.0 / 8.0 * hhmin * qual->mu_min / (xnorm + hhmin);
	//k *= 4.0;
	//k = 0.01;
	//	k = 1.0e-3;
        //k = 1.0e-1;
	//k = 1.0 / 2.0 * hhmin * qual->mu_min / (xnorm + hhmin);
	//k = 1.0 * 1.0 / 2.0 * hhmin * qual->mu_min / (xnorm + 1.0e-7);


	cout << "smoother_incr: " << incr << endl;
	cout << "t: " << t << endl;
	cout << "T: " << T << endl;
	cout << "num_steps: " << num_steps << endl;
	
	num_steps++;
	korig = k;

	//message("ElasticSmoother timer step: %g", MPI::stopTimer(timer0));


	if(num_steps >= max_steps)
	{
	  return false;
	}
	return true;

      }

      static void deform(Mesh& mesh, real k, Function& u);
      void deform_x(Function& X);

      Form* aS;
      Form* LS;
      Form* aP;
      Form* LP;

      Function& U;
      Function& B;
      Function& B0;
      TimeStep_Function& kf;
      Function& icv;
      Function X;
      Function X0;
      Vector B0x, Bx, Btmp;
      Matrix BM;
      Vector Bm;
      Vector X0x, Xx;
      Vector icell_vol;
      
      Vector xx0, xxtmp;
      
      real* vxvalues;

      real hmin;
      real kk;
      real uhmin0;
      MeshFunction<real>& h0;

      MeshQuality* qual;

      int num_steps;
      int max_steps;
      int min_steps;

//       real timer0;
//       real timer1;
//       real timer2;
    };
  };
//-----------------------------------------------------------------------------
}}

#endif
