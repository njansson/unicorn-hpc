#ifndef __ERROR_ESTIMATE_H
#define __ERROR_ESTIMATE_H

#include <dolfin.h>

namespace dolfin {
  namespace unicorn {
    class ErrorEstimate
    {
    public:
      
      /// Constructor
      ErrorEstimate(Mesh& mesh, Form* Lres_1, Form* Lgradphi = 0);

      /// Constructor
      ErrorEstimate(Mesh& mesh, Form* Lres_1, Form* Lres_2, 
		    Form* Lres_3, Form* Lgradphi = 0);

      /// Destructor
      ~ErrorEstimate();
      
      /// Initialize the class
      void init(Mesh& mesh, Form* Lres_1, Form* Lres_2, Form* Lres_3, 
		Form* Lgradphi, MeshFunction<real>& Rf, MeshFunction<real>& eif) ;

      // Compute error (norm estimate)
      void ComputeError(real& error);
      
      // Compute error indicator
      void ComputeErrorIndicator(real t, real k, real T, real w = 1.0);
      
      // Compute largest indicators
      void ComputeLargestIndicators(std::vector<int>& cells, real percentage);

      // Refine based on indicators
      void AdaptiveRefinement(real percentage);
      
      // Comparison operator for index/value pairs
      struct less_pair : public std::binary_function<std::pair<int, real>,
	std::pair<int, real>, bool>
      {
	bool operator()(std::pair<int, real> x, std::pair<int, real> y)
	{
	  return x.second < y.second;
	}
      };


      Mesh& mesh;
      
      // For full conservation law three are residuals
      Form* Lres_1;
      Form* Lres_2;
      Form* Lres_3;
      
      Form* Lgradphi;
      Function gradphi;
      Function res;

      // Error indicator
      Vector e_indx;

      // Residual
      Vector res_1x, res_2x, res_3x;

      // Gradient of phi
      Vector gradphix;
      Assembler assembler;

      // Norm of residual
      MeshFunction<real> Rf;
      // Error indicator
      MeshFunction<real> eif;      

    private:

      // Compute largest indicators (eind based)
      void ComputeLargestIndicators_eind(std::vector<int>& cells, real percentage);

      // Compute largest indicators (cell based)
      void ComputeLargestIndicators_cell(std::vector<int>& cells, real percentage);
      
      void merge(real *a,real *b,real *res,int an,int bn);
    };
    
  }}

#endif
