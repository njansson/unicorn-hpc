#ifndef __ELASTIC_SMOOTHER_H
#define __ELASTIC_SMOOTHER_H

#include <dolfin.h>
#include <unicorn/EquiAffineMap.h>
#include <unicorn/MeshQuality.h>
#include <unicorn/TimeDependentPDE.h>
#include <unicorn/Project.h>

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
		MeshFunction<real>& h0,
		Function& W, real k);

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

	//real p = 2.0;
	real p = 2.0;

	real val = 1.0;

	real qq = mqual.cellQuality(c);
	real qual = 1.0 / pow(qq, p);
	//real qual = q.get(cell().index());
	//real scale = 1.0 * maxqual / minqual;
	real scale = 1.0;
	val = qual / scale;
 	//val = 1.0;

	values[0] = val;
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
      real scale = h;
      //real scale = pow(cell.volume(), 1.0 / cell.dim()) / cell.diameter() * h;



      real Finv[3*3];
      //      uBlasDenseMatrix Finv(map.C);
      //      Finv *= scale;
      for (uint ii = 0; ii < 3; ii++)
	for (uint jj = 0; jj < 3; jj++)
	  Finv[RM(ii,jj,3)] = map.C[RM(ii,jj,3)] * scale;
      

      
      real B[3*3];
      //      uBlasDenseMatrix B;
      //      B.mat() = ublas::prod(ublas::trans(Finv.mat()), Finv.mat());
      for (uint ii = 0; ii < 3; ii++) 
	for(uint jj = 0; jj < 3; jj++)
	  for (uint r =0; r < 3; r++)
	    B[RM(ii,jj,3)] += (Finv[RM(r,ii,3)] * Finv[RM(r,jj,3)]);

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
	icv.init(mesh, icell_vol, *LS, 3);
	

	//     B.init(mesh, Bx, *L, 2);
	//     B0.init(mesh, B0x, *L, 2);
	
	// X is compatible with U
	X.init(mesh, Xx, *a, 0);
	//X0x.init(X.vector().size());
	X0.init(mesh, X0x, *a, 0);

	int d = mesh.topology().dim();
	
	vxvalues = new real[d * mesh.numVertices()];

	qual = new MeshQuality(mesh);

	xx0.init(X0.vector().local_size());
	xxtmp.init(X0.vector().local_size());
	//xx0 = 0.0;

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

	assembler->assemble(Btmp, *LS);
	dolfin_set("output destination", "terminal");
		
	B.vector() = Btmp;
	
	B.vector() *= k;
	B.vector()+= B0.vector();
	message("computeB took %g seconds",toc());
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
      }
      
      void preparestep()
      {
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
// 	cout << "Ubefore: " << endl;
// 	U.vector().vec().disp();
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
	//deform(mesh(), X);

//  	cout << "X0 1: " << endl;
//  	X0.vector().disp();

	qual->meshQuality();

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
	
	dolfin_set("output destination", "silent");
	assembler->assemble(dotx, L(), reset_tensor);
	//     for (uint i = 0; i < bc_mom.size(); i++)
	//       bc_mom[i]->apply(M, dotx, a());
	dolfin_set("output destination", "terminal");
// 	set("output destination", "silent");
// 	assembler.assemble(dotx, L());
// 	for (uint i = 0; i < bc().size(); i++)
// 	  bc()[i]->apply(M, dotx, a());
// 	set("output destination", "terminal");
	
	//cout << "dtU: " << dotx.norm(linf) << endl;
      }
      
      void save(Function& U, real t)
      {
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
      }

      void u0(Vector& x)
	//void u0(uBlasVector& x)
      {
	cout << "ElasticityPDE::u0" << endl;

	int d = mesh().topology().dim();
	int N = mesh().numVertices();
	int M = mesh().numCells();

	//	real* Xarr = X.vector().vec().array();
	real* Xarr = new real[X.vector().local_size()];

	for (VertexIterator n(mesh()); !n.end(); ++n)
	{
	  for(int i = 0; i < d; i++)
	  {
	    Xarr[i * N + n->index()] = n->x()[i];
	  }
	}
	
	X.vector().set(Xarr);
	X.vector().apply();
	delete[] Xarr;

	U.vector().zero();

	//	real* Barr = B.vector().vec().array();
	real* Barr = new real[B.vector().local_size()];
	
	for (CellIterator c(mesh()); !c.end(); ++c)
	{
	  for(int i = 0; i < d * d; i++)
	  {
	    real val = initial_B(*c, h0, i);
	    
	    Barr[i * M + c->index()] = val;
	  }
	}

	B.vector().set(Barr);
	B.vector().apply();
// 	// Compute initial value
// 	Function Pf = pproject(*aP, *LP, mesh());

// 	U.vector().vec().copy(Pf.vector().vec(), 0, 0, U.vector().size());

// 	cout << "u0: " << endl;
// 	Pf.vector().disp();
	
// 	// Copy into components
// 	X.vector().vec().copy(Pf.vector().vec(), 0, 0, X.vector().size());
// 	U.vector().vec().copy(Pf.vector().vec(), 0, U.vector().size(),
// 			      U.vector().size());
// 	B.vector().vec().copy(Pf.vector().vec(), 0, 2*U.vector().size(),
// 			      B.vector().size());

	
	
	//FIXME: shift() seems to be needed for cG(1) but not dG(0), why?
	shift();
	
//  	cout << "X0: " << endl;
//  	X0.vector().disp();
// 	cout << "x0: " << endl;
// 	x.disp();
// 	cout << "B0: " << endl;
// 	B.vector().disp();
      }

      //bool update(const uBlasVector& u, real t, bool end)
      bool update(real t, bool end)
      {
 	std::cout << "ElasticSmoother::update:" << std::endl;
// 	cout << "U:" << endl;
// 	U.vector().disp();
// 	cout << "B:" << endl;
// 	B.vector().disp();

	deform(mesh(), k, U);

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
	cout << "dotx norm: " << norm << endl;

	real xnorm = U.vector().norm(linf);
	std::cout << "x norm: " << xnorm << std::endl;
	//u.disp();

	xxtmp = xx0;
	xxtmp -= U.vector();
	real xinc = xxtmp.norm(linf) / xnorm;

	static real xn0;
	cout << "x inc: " << xinc << endl;

	real xn = xnorm / qual->h_min;
	cout << "xn: " << xn << endl;

	cout << "hmin: " << qual->h_min << endl;

	real xmax = std::max(x->norm(linf), 1.0e-6);

	//	real geom_scale = qual->bbox_min.distance(qual->bbox_max);

	real uhmin = 1.0e9;
	real vumax = 0.0;
	real hhmin = 1.0e9;

	//	real* z = U.vector().vec().array();
	real* z = new real[U.vector().local_size()];
	for(CellIterator ci(mesh()); !ci.end(); ++ci)
	{
	  for(EdgeIterator ei(*ci); !ei.end(); ++ei)
	  {
	    Vertex v0(mesh(), ei->entities(0)[0]);
	    Vertex v1(mesh(), ei->entities(0)[1]);
	    
	    real vrel = 0.0;
	    for(int i = 0; i < d; i++)
 	    {
	      vrel += sqr(z[i * N + v0.index()] - z[i * N + v1.index()]);
	    }
	    vrel = sqrt(vrel);

	    real h = pow(ci->volume(), 1.0 / ci->dim());

 	    real qi = qual->cellQuality(*ci);
 	    real uhi = qi * h / vrel;

 	    uhmin = std::min(uhmin, uhi);
 	    vumax = std::max(vumax, vrel);
 	    hhmin = std::min(hhmin, h);

	    //wmax = std::max(wmax, vrel);
	  }

// 	  for (VertexIterator vi(*ci); !vi.end(); ++vi)
// 	  {
// 	    Vertex& v = *vi;
	    
// 	    real vu = 0.0;
// 	    for ( int i = 0; i < d; i++ )
// 	    {
// 	      vu += sqr(u[i*N + v.index()]);
// 	    }
// 	    vu = sqrt(vu);
	    
// 	    real qi = qual->cellQuality(*ci);
// 	    real uhi = qi * ci->diameter() / vu;

// 	    uhmin = std::min(uhmin, uhi);
// 	    vumax = std::max(vumax, vu);
// 	    hhmin = std::min(hhmin, ci->diameter());
// 	  }
	}

	U.vector().set(z);
	U.vector().apply();

	delete[] z;

	if(t == 0.0)
	  uhmin0 = uhmin;

	real uh = std::min(uhmin, uhmin0);
	
	cout << "uhmin: " << uhmin << endl;
	cout << "vumax: " << vumax << endl;
	cout << "hhmin: " << hhmin << endl;

	//real kcfl = 0.1 * qual->h_min * qual->mu_min / xmax;
	real kcfl = 0.01 * uh;
	real kcons = 0.1 * qual->h_min *qual->mu_min / std::max(xmax, 1.0);
	//kcfl = 2.0 * kcfl * kk0 / (kcfl + kk0);
	//k = (k + kk0) / 2.0;

	uhmin0 = uhmin;


	// Compute k according to k < h / |U|
// 	if(true || t == 0.0)
// 	  k = kk;
// 	else
// 	  k = 0.25 * qual->h_min * qual->mu_min / xmin;

	if(t == 0.0)
	  k = std::min(kk, kcfl);
	else
	  k = kcfl;

	//real ramp = 4.0 * (1.0 + t / (T * 10.0));
	//real ramp = pow(1.0 + t / (T / 1000.0), 5.0);
	//real ramp = pow(1.0 + 1.0e2 * 10.0 * t, 3.0);
	static real ramp;

	if(t == 0.0)
	  ramp = 2.0;
	else
	  ramp *= 2.0;

	if(t == 0.0)
	  xn0 = 1.0e3;
	else if(xn > xn0)
	{
	  cout << "diverging: " << xn << " " << xn0 << endl;
	  ramp /= 10.0;
	}

	//k *= ramp;

// 	if(fabs(k - 1.0) < 0.2)
// 	  k = 0.8;

	//k = std::min(k, kcfl);
	k = std::min(k, 1.0e6);
	k = std::min(k, ramp * kcons);
	if(fabs(k - 1.0) < 0.1)
	  k = 0.5;

	//k = std::min(k, 1.0e-5);
	
        k = 1.0 / 2.0 * hhmin * qual->mu_min / (xnorm + hhmin);
	//	k = 1.0e-3;
        //k = 1.0e-1;
	//k = 1.0 / 2.0 * hhmin * qual->mu_min / (xnorm + hhmin);
	//k = 1.0 * 1.0 / 2.0 * hhmin * qual->mu_min / (xnorm + 1.0e-7);


	cout << "kcfl: " << kcfl << endl;
	cout << "t: " << t << endl;
	cout << "T: " << T << endl;
	cout << "num_steps: " << num_steps << endl;
	
	num_steps++;
	xx0 = U.vector();

	xn0 = xn;

	korig = k;

	//if(num_steps > min_steps &&
	//   (xinc < sstate_tol || num_steps > max_steps))
	if((num_steps > min_steps &&
	    (xn < sstate_tol || num_steps > max_steps)))
	{
	  return false;
	}
	return true;

      }

      static void deform(Mesh& mesh, real k, Function& u);

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
