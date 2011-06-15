// Copyright (C) 2008 Murtazo Nazarov
// Licensed under the GNU GPL Version 2.
//
// Modified by Niclas Jansson 2009-2011.
//
#include <dolfin.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>


#include "unicorn/UniParameters.h"

using namespace dolfin;
using namespace dolfin::unicorn;

//-----------------------------------------------------------------------------
void UniParameters::parse_parameters(std::string fname)
{

  std::map<std::string, Type> parameters;
  
  parameters["alpha"] = _REAL_;
  parameters["beta"] = _REAL_;
  parameters["nu"] = _REAL_;
  parameters["Ma"] = _REAL_;
  parameters["n_samples"] = _INT_;
  parameters["adapt_iter"] = _INT_;
  parameters["adapt_tol"] = _REAL_;
  parameters["adapt_percentage"] = _REAL_;
  parameters["adapt_algorithm"] = _STR_;
  parameters["adapt_project"] = _BOOL_;
  parameters["adapt_type"] = _STR_;
  parameters["T"] = _REAL_;
  parameters["Ubar"] = _REAL_;
  parameters["dual_T"] = _REAL_;
  parameters["krylov_method"] = _STR_;
  parameters["krylov_pc"] = _STR_;
  parameters["krylov_pc_keep"] = _BOOL_;
  parameters["output_format"] = _STR_;

  std::ifstream param_file;

  std::string str;

  param_file.open(fname.c_str());

  if (!param_file.good())
    error("The parameters file %s cannot be opened", fname.c_str());

  const char *sep = " -\n\t";
  char *token;

  while (param_file.good()) 
  {
    getline(param_file, str);
    std::string str_value;
    uint state = 0;
    for(token = strtok(const_cast<char *>(str.c_str()), sep); 
	token; token = strtok(NULL, sep), (++state %= 2))
    {
      if ( state == 0)
      {
	str_value = token;
	continue;
      }
      std::map<std::string, Type>::iterator it = parameters.find(token);

      if ( it != parameters.end())
      {
	if( it->second  == _INT_ )
	{
	  int value;
	  if (parse_numeric(&value, str_value))
	    dolfin_add(token, value);
	  else
	    error("Failed to parse numeric type");
	}
	else if( it->second  == _REAL_ )
	{
	  real value;
	  if (parse_numeric(&value, str_value))
	    dolfin_add(token, value);
	  else
	    error("Failed to parse numeric type");
	}
	else if( it->second  == _STR_ )
	{
	  dolfin_add(token, str_value);
	}
	else if( it->second == _BOOL_ )
	{
	  int value;
	  if (parse_numeric(&value, str_value))
	  {
	    if (value)
	      dolfin_add(token, true);
	    else
	      dolfin_add(token, false);
	  }
	  else
	    error("Failed to parse numeric type");
	}
      }
    }        
  }
  param_file.close();
}
//-----------------------------------------------------------------------------

