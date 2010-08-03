#ifndef __ELASTIC_SMOOTHER_H
#define __ELASTIC_SMOOTHER_H

#include <dolfin.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <unicorn/EquiAffineMap.h>
#include <unicorn/MeshQuality.h>
#include <unicorn/TimeDependentPDE.h>
#include <unicorn/Project.h>
#include <dolfin/fem/UFC.h>

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
		MeshFunction<bool>& masked_cells,
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
    Quality(Mesh& mesh, MeshFunction<real>& h0) : mqual(mesh), Function(mesh)
      {
      }
    
      void eval(real* values, const real* x) const
      {
	//FIXME: should we have to cast?
	Cell& c = const_cast<Cell&>(cell());

	//real p = 2.0;
	real p = 1.0;

	real val = 1.0;

	real qq = mqual.cellQuality(c);
	real qual = 1.0 / pow(qq, p);
	real scale = 1.0;
	val = qual / scale;

	values[0] = val;
      }

      MeshQuality mqual;
    };

    static real initial_B(Cell& cell, MeshFunction<real>& h0, int i)
    {
      EquiAffineMap map;
      map.update(cell);
      
      real h = h0.get(cell.index());
      real scale = h / sqrt(cell.dim());

      uBlasDenseMatrix Finv(3, 3);
      for (uint ii = 0; ii < 3; ii++)
	for (uint jj = 0; jj < 3; jj++)
	{
 	  Finv.mat()(ii, jj) = map.C[RM(ii,jj,3)] * scale;
	}
      uBlasDenseMatrix B(3, 3);
      B.mat() = ublas::prod(ublas::trans(Finv.mat()), Finv.mat());

      int d = cell.dim();

      int ii = i % d;
      int jj = i / d;

      return B.mat()(ii, jj);
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
	//	uBlasDenseMatrix Finv(map.C);
	//	Finv *= scale;
	for (uint ii = 0; ii < 3; ii++)
	  for (uint jj = 0; jj < 3; jj++)
	    Finv[RM(ii,jj,3)] = map.C[RM(ii,jj,3)] * scale;
      

	real B[3*3];
//        uBlasDenseMatrix B;
//	B.mat() = ublas::prod(ublas::trans(Finv.mat()), Finv.mat());
	for (uint ii = 0; ii < 3; ii++) 
	  for(uint jj = 0; jj < 3; jj++)
	    for (uint r =0; r < 3; r++)
	      B[RM(ii,jj,3)] += (Finv[RM(r,ii,3)] * Finv[RM(r,jj,3)]);



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
    MyBC(Mesh& mesh) : Function(mesh)
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
      TimeDependentPDE(mesh, bc, T),
	aS(aS), LS(LS),
	U(U), B(B), B0(B0), kf(kf), icv(icv),
	hmin(0.0), kk(kk), h0(h0),
	qual(0), num_steps(0), max_steps(dolfin_get("Smoother max time steps")),
	min_steps(4)
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
	
	dolfin_set("output destination", "silent");
	//assembler.assemble(BM, *aS);
	//BM.mat().lump(Bm.vec());

	ComputeInvCellVolume(icell_vol);

	Btmp.zero();
	assembler->assemble(Btmp, *LS);
	dolfin_set("output destination", "terminal");
		
	B.vector() = Btmp;
	
	B.vector() *= k;
	B.vector()+= B0.vector();
	if(dolfin::MPI::processNumber() == 0)
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
	X0.vector().apply();
	B0.vector().apply();
      }
      
      void preparestep()
      {
	shift();
      }
      
      void prepareiteration()
      {
	computeX();
	computeB();
	deform_x(X);

	qual->meshQuality();

	cout << "U: " << x->norm(linf) << endl;

	reassemble = true;
      }
      
      void fu(const Vector& x, Vector& dotx, real T)
      {
	dolfin_set("output destination", "silent");
	assembler->assemble(dotx, L(), reset_tensor);
	dolfin_set("output destination", "terminal");
      }
      
      void save(Function& U, real t)
      {
	if(dolfin::MPI::processNumber() == 0)
	  cout << "Saving at t: " << t << endl;

	if(t == 0.0)
	{
	  solutionfile << U;
	}
	
	if(true)
	  //while(lastsample + sampleperiod < t)
	//if(lastsample + sampleperiod < t)
	{
	  if(dolfin::MPI::processNumber() == 0)
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
      
      
	  if(dolfin::MPI::processNumber() == 0)
	    qual->disp();
	}
      }

      void u0(GenericVector& x)
      {
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

      }

      void revert()
      {
	B.vector() = B0.vector();
	B.vector().apply();
      }

      bool update(real t, bool end)
      {
	if(dolfin::MPI::processNumber() == 0)
	  std::cout << "ElasticSmoother::update:" << std::endl;

	qual->meshQuality();
	//qual.disp();
	uint inverted_cell;
	
	if(qual->isInverted(inverted_cell)) {
	  cout << "inverted cell: " << inverted_cell << endl;
	}

	int d = mesh().topology().dim();
	int N = mesh().numVertices();

	real sstate_tol = 1.0e-9;

	real norm = dotx->norm(linf);

	real xnorm = U.vector().norm(linf);

	real hhmin = 1.0e9;

	for(CellIterator ci(mesh()); !ci.end(); ++ci)
	{
	  real h = pow(ci->volume(), 1.0 / ci->dim());
	  real qi = qual->cellQuality(*ci);

	  hhmin = std::min(hhmin, h);
	}
	hhmin = MeshQuality::reduce_min_scalar(hhmin);


        k = 1.0 / 2.0 * hhmin * qual->mu_min / (xnorm + hhmin);

	if(dolfin::MPI::processNumber() == 0)
	{
	  cout << "dotx norm: " << norm << endl;
	  cout << "x norm: " << xnorm << endl;
	  cout << "hhmin: " << hhmin << endl;
	  
	  cout << "t: " << t << endl;
	  cout << "T: " << T << endl;
	  cout << "num_steps: " << num_steps << endl;
	}

	
	num_steps++;
	korig = k;

	if((num_steps > min_steps &&
	    (num_steps > max_steps)))
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
    };
  };
//-----------------------------------------------------------------------------
}}

#endif
