// Copyright (C) 2008 Murtazo Nazarov
// Licensed under the GNU GPL Version 2.
//
// Modified by Niclas Jansson 2009-2010.
//
#include <dolfin.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <string>
#include <boost/tokenizer.hpp>


#include "unicorn/UniParameters.h"

typedef boost::tokenizer<boost::char_separator<char> > tokenizer;   
using namespace dolfin;

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
  parameters["adapt_type"] = _STR_;
  parameters["T"] = _REAL_;
  parameters["Ubar"] = _REAL_;
  parameters["dual_T"] = _REAL_;
  parameters["krylov_method"] = _STR_;
  parameters["krylov_pc"] = _STR_;

  std::ifstream param_file;

  std::string str;

  param_file.open(fname.c_str());

  if (!param_file.good())
    error("The parameters file %s cannot be opened", fname.c_str());

  boost::char_separator<char> sep(" -\n\t");

  while (param_file.good()) 
  {
    getline(param_file, str);
    tokenizer tok(str,sep);
    std::string str_value;
    uint state = 0;
    for(tokenizer::iterator beg=tok.begin(); beg!=tok.end(); beg++,(++state %= 2))
    {
      if ( state == 0)
      {
	str_value = *beg;
	continue;
      }
      std::map<std::string, Type>::iterator it = parameters.find(*beg);

      if ( it != parameters.end())
      {
	if( it->second  == _INT_ )
	{
	  int value;
	  if (parse_numeric(&value, str_value))
	    dolfin_add(*beg, value);
	  else
	    error("Failed to parse numeric type");
	}
	else if( it->second  == _REAL_ )
	{
	  real value;
	  if (parse_numeric(&value, str_value))
	    dolfin_add(*beg, value);
	  else
	    error("Failed to parse numeric type");
	}
	else if( it->second  == _STR_ )
	{
	  dolfin_add(*beg, str_value);
	}
      }
    }        
  }
  param_file.close();
}
//-----------------------------------------------------------------------------

