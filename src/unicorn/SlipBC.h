// Copyright (C) 2007 Murtazo Nazarov
// Licensed under the GNU LGPL Version 2.1.
//
// Existing code for Dirichlet BC is used
//
// Modified by Niclas Jansson, 2008-2009.
//
// First added:  2007-05-01
// Last changed: 2010-05-09

#ifndef __SLIPBC_H
#define __SLIPBC_H

#include <dolfin.h>
#include "NodeNormal.h"



namespace dolfin
{

  class DofMap;
  class Function;
  class Mesh;
  class SubDomain;
  class Form;
  class GenericMatrix;
  class GenericVector;
  
  class SlipBC : public BoundaryCondition
  {
  public:
    
    /// Create boundary condition for sub domain
    SlipBC(Mesh& mesh, SubDomain& sub_domain);

    /// Create boundary condition for sub domain
    SlipBC(Mesh& mesh, SubDomain& sub_domain, NodeNormal& Node_normal);
    
    /// Create boundary condition for sub domain specified by index
    SlipBC(MeshFunction<uint>& sub_domains, uint sub_domain);
    
    /// Create sub system boundary condition for sub domain
    SlipBC(Mesh& mesh, SubDomain& sub_domain, const SubSystem& sub_system);
    
    /// Create sub system boundary condition for sub domain specified by index
    SlipBC(MeshFunction<uint>& sub_domains, uint sub_domain, const SubSystem& sub_system);

    /// Destructor
    ~SlipBC();

    /// Apply boundary condition to linear system
    void apply(GenericMatrix& A, GenericVector& b, const Form& form);
    
    void apply(GenericMatrix& A, GenericVector& b, const GenericVector& x, 
	       const DofMap& dof_map, const ufc::form& form); 

    void apply(GenericMatrix& A, GenericVector& b, const GenericVector& x, 
	       const Form& form);
    



  private:

    void applySlipBC(Matrix& A, Matrix& As, Vector&, Mesh& mesh, uint node, Array<uint>& nodes);
    
    // Do: A(row,col) = value   using setblock not setvalue
    inline void Aset(Matrix& A, uint row, uint col, real value) { A.set(&value, 1, &row, 1, &col);};

    // Do: b(row) = value   using setblock not setvalue
    inline void bset(Vector& b, uint row, real value) { b.set(&value, 1, &row); };


    // Initialize sub domain markers    
    void init(SubDomain& sub_domain);


    void apply(GenericMatrix& A, GenericVector& b, const DofMap& dof_map, 
	       const ufc::form& form);


    void apply(GenericMatrix& A, GenericVector& b, const DofMap& dof_map, 
	       const Form& form);

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

    // User defined sub domain
    SubDomain* user_sub_domain;

    // Node normal and tangents
    NodeNormal node_normal;

    int nzm;

    Matrix* As;

    int N_local;
    int N_offset;
    std::set<uint> off_proc_rows;

    real *row_block;
    real *zero_block;
    uint *a1_indices_array;

    BoundaryMesh* boundary;
    MeshFunction<uint> *cell_map, *vertex_map;
  };
}

#endif
