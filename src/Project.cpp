// Copyright (C) 2006 Johan Jansson.
// Licensed under the GNU GPL Version 2.
//
// Modified by Anders Logg 2006.
//
// First added:  2006
// Last changed: 2006-05-04

#include <dolfin.h>
#include <unicorn/Project.h>

using namespace dolfin;
using namespace dolfin::unicorn;
//-----------------------------------------------------------------------------
Function dolfin::unicorn::pproject(Form& aP, Form& LP, Mesh& mesh)
{
  Matrix M;
  Vector b;
  Vector *x = new Vector();

  Function Pf(mesh, *x, aP, 1);
  
  tic();
  dolfin_set("output destination", "silent");
  Assembler assembler(mesh);
  assembler.assemble(M, aP);
  assembler.assemble(b, LP);
  dolfin_set("output destination", "terminal");
  message("Project assemble took %g seconds",toc());
    
  tic();
  KrylovSolver solver(bicgstab, amg);
  solver.solve(M, *x, b);
  message("Project linear solve took %g seconds",toc());

  cout << "pproject:" << endl;
  Pf.vector().disp();
  
  return Pf;
}
//-----------------------------------------------------------------------------

