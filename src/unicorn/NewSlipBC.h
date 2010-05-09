// Copyright (C) 2007 Murtazo Nazarov
// Licensed under the GNU LGPL Version 2.1.
//
// Existing code for Dirichlet BC is used
//
// First added:  2007-05-01
// Last changed: 2008-04-01

#ifndef __NEWSLIPBC_H
#define __NEWSLIPBC_H

#include <dolfin/common/constants.h>
#include <dolfin/fem/SubSystem.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/SubDomain.h>
#include <dolfin/fem/UFCMesh.h>
#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/BoundaryCondition.h>
#include <dolfin/la/Matrix.h>

#include <unicorn/NodeNormal.h>

namespace dolfin
{

  class DofMap;
  class Function;
  class Mesh;
  class SubDomain;
  class Form;
  class GenericMatrix;
  class GenericVector;
  class NormalNode;
  
  class NewSlipBC : public BoundaryCondition
  {
  public:
    
    /// Create boundary condition for sub domain
    NewSlipBC(Mesh& mesh, SubDomain& sub_domain);
    
    /// Create boundary condition for sub domain specified by index
    NewSlipBC(MeshFunction<uint>& sub_domains, uint sub_domain);
    
    /// Create sub system boundary condition for sub domain
    NewSlipBC(Mesh& mesh, SubDomain& sub_domain, const SubSystem& sub_system);
    
    /// Create sub system boundary condition for sub domain specified by index
    NewSlipBC(MeshFunction<uint>& sub_domains, uint sub_domain, const SubSystem& sub_system);

    /// Destructor
    ~NewSlipBC();
    
    /// Apply boundary condition to linear system
    void apply(GenericMatrix& A, GenericVector& b, const Form& form);
    
    void apply(GenericMatrix& A, GenericVector& b, const GenericVector& x, 
	       const DofMap& dof_map, const ufc::form& form); 
    
    void apply(GenericMatrix& A, GenericVector& b, const GenericVector& x, 
	       const Form& form);
    
    void applyNewSlipBC(Matrix& A, Matrix& As, Vector& b, Mesh& mesh, 
			uint node, Array<uint>& nodes);
    
    // Do: A(row,col) = value   using setblock not setvalue
    void Aset(Matrix& A, int row, int col, real value);

    // Do: b(row) = value   using setblock not setvalue
    void bset(Vector& b, int row, real value);

    std::vector< std::vector<int> > permutations;
    

    std::vector< std::vector<int> > permutations;

  private:

    // Initialize sub domain markers    
    void init(SubDomain& sub_domain);


    void apply(GenericMatrix& A, GenericVector& b, const DofMap& dof_map, 
	       const ufc::form& form);
    
    void apply(GenericMatrix& A, GenericVector& b, const DofMap& dof_map, 
	       const Form& form);



    // Do: A(row,col) = value   using setblock not setvalue
    inline void Aset(Matrix& A, uint row, uint col, real value) { A.set(&value, 1, &row, 1, &col);};

    // Do: b(row) = value   using setblock not setvalue
    inline void bset(Vector& b, uint row, real value) { b.set(&value, 1, &row); };

    // The mesh
    Mesh& mesh;

    // Sub domain markers (if any)
    MeshFunction<uint>* sub_domains;

    // The sub domain
    uint sub_domain;

    // True if sub domain markers are created locally
    bool sub_domains_local;

    // Sub system
    SubSystem sub_system;

    // Node normal and tangents
    NodeNormal node_normal;

    // User defined sub domain
    SubDomain* user_sub_domain;



    int nzm;

    Matrix* As;


    int N_local;
    int N_offset;
    std::set<uint> off_proc_rows;

    real *row_block;
    real *zero_block;
    uint *a1_indices_array;

  };
}

#endif
