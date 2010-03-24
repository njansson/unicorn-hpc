// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU GPL Version 2.
//
// First added:  2005-11-26
// Last changed: 2006-12-12

#ifndef __SPACE_TIME_FUNCTION_H
#define __SPACE_TIME_FUNCTION_H

#include <ufc.h>
#include <dolfin/function/GenericFunction.h>

namespace dolfin { namespace unicorn
{
  //  class dolfin::Function;

  class SpaceTimeFunction
  {
  public:

    /// Create space-time function
    SpaceTimeFunction(Mesh& mesh, Function& Ut);

    /// Copy constructor
    SpaceTimeFunction(const SpaceTimeFunction& f);

    /// Destructor
    ~SpaceTimeFunction();

    /// Evaluate function at time t, giving result in Ut
    void eval(real t);

    // Add a space function at time t
    void addPoint(std::string Uname, real t);

    void util_addFiles(std::vector<std::string> filenames, real T);
    void util_fileList(std::string basename, int N,
		       std::vector<std::string>& filenames);

    /// Return mesh associated with function
    Mesh& mesh();

    /// Return interpolant function
    Function& evaluant();

  protected:

    // Pointer to mesh associated with function (null if none)
    Mesh* _mesh;

    // Pointer to evaluant function
    Function* _function;

    // Space functions defining the current time interval (cache)
    Function U0;
    Function U1;
    //    Function* dtU0;
    //    Function* dtU1;

    std::map<real, std::string> U_files;
  };

}}

#endif
