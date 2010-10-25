// Copyright (C) 2008 Murtazo Nazarov
// Licensed under the GNU GPL Version 2.
//
// Modified by Niclas Jansson 2009-2010.
//
#ifndef __UNIPARAMETERS_H
#define __UNIPARAMETERS_H

#include <string>
#include <sstream>
#include <dolfin.h>


namespace dolfin
{
  namespace unicorn 
  {
    class UniParameters 
    {
    public:
      
      static void parse_parameters(std::string);    
      
    private:
      
      enum Type { _INT_, _REAL_, _STR_, _BOOL_};
      
      template<typename T>
      static bool parse_numeric(T* i, std::string s)
      {
	std::istringstream ss(s);
	if(ss>>*i)
	  return true;
	else
	  return false;
      }        
    };
  }
}

#endif
