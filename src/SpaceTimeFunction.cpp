// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU GPL Version 2.
//
// First added:  2005-11-26
// Last changed: 2006-12-12
//
// Note: this breaks the standard envelope-letter idiom slightly,
// since we call the envelope class from one of the letter classes.

#include <sstream>
#include <dolfin/la/Vector.h>
#include <dolfin/function/Function.h>
#include <dolfin/io/File.h>
#include "unicorn/SpaceTimeFunction.h"
#include <dolfin/main/MPI.h>

using namespace dolfin;
using namespace dolfin::unicorn;

//-----------------------------------------------------------------------------
SpaceTimeFunction::SpaceTimeFunction(Mesh& mesh, Function& Ut)
  : _mesh(&mesh), _function(&Ut),
    U0(Ut), U1(Ut)
{
}
//-----------------------------------------------------------------------------
SpaceTimeFunction::SpaceTimeFunction(const SpaceTimeFunction& f)
{
}
//-----------------------------------------------------------------------------
SpaceTimeFunction::~SpaceTimeFunction()
{
}
//-----------------------------------------------------------------------------
void SpaceTimeFunction::eval(real t)
{
  // FIXME: Implement a cache

  // Find last element
  //  real T = (*(U_files.rbegin())).first;
//   cout << "eval t: " << t << endl;
//   cout << "T: " << T << endl;

//  dolfin_assert(t <= T);

//   std::map<real, std::string>::iterator it1 = U_files.upper_bound(t);
//   std::map<real, std::string>::iterator it0 = it1;
//   --it0;

  std::map<real, std::string>::iterator it1;
  std::map<real, std::string>::iterator it0;

  // Find element in U_files so that element < t
  it1 = U_files.upper_bound(t);

  // If t == T, we need to step back one
  if(it1 == U_files.end())
    --it1;

  it0 = it1;
  --it0;
  
  real t0 = (*it0).first;
  real t1 = (*it1).first;

  cout << "t0: " << t0 << endl;
  cout << "t1: " << t1 << endl;

  std::string name0 = (*it0).second;
  std::string name1 = (*it1).second;

  cout << "name0: " << name0 << endl;
  cout << "name1: " << name1 << endl;

  File file0(name0);
  File file1(name1);


  file0 >> U0.vector();
  file1 >> U1.vector();


  // Compute weights (linear Lagrange interpolation)
  real w0 = (t1 - t) / (t1 - t0);
  real w1 = (t - t0) / (t1 - t0);

  // Compute interpolated value
  evaluant().vector() = 0.0;
  evaluant().vector().axpy(w0, U0.vector());
  evaluant().vector().axpy(w1, U1.vector());
}
//-----------------------------------------------------------------------------
void SpaceTimeFunction::addPoint(std::string Uname, real t)
{
  U_files[t] = Uname;
}
//-----------------------------------------------------------------------------
void SpaceTimeFunction::util_addFiles(std::vector<std::string> filenames,
				      real T)
{
  //FIXME: For now we assume a fixed time step

  int counter = 0;
  int num_files = filenames.size();

  for(std::vector<std::string>::iterator it = filenames.begin();
      it != filenames.end(); ++it)
  {
    std::string filename = *it;

    real t = T * real(counter) / real(num_files - 1);

    addPoint(filename, t);

    counter++;
  }
}
//-----------------------------------------------------------------------------
void SpaceTimeFunction::util_fileList(std::string basename, int N,
				      std::vector<std::string>& filenames)
{
  filenames.clear();

  for(int i = 0; i < N - 1; i++)
  {
    std::stringstream filename, number;
    number.fill('0');
    number.width(6);

    number << i;

    filename << basename;
    filename << number.str();
    filename <<  "_" << MPI::processNumber() << ".bin";
    filename << std::ends;

    filenames.push_back(filename.str());
  }
}
//-----------------------------------------------------------------------------
Mesh& SpaceTimeFunction::mesh()
{
  dolfin_assert(_mesh);
  return *_mesh;
}
//-----------------------------------------------------------------------------
Function& SpaceTimeFunction::evaluant()
{
  dolfin_assert(_function);
  return *_function;
}
//-----------------------------------------------------------------------------
