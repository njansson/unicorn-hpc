// Copyright (C) 2007 Murtazo Nazarov
// Licensed under the GNU LGPL Version 2.1.
//
// Existing code for Dirichlet BC is used
//
// Modified by Niclas Jansson, 2008-2010.
// 
// First added:  2007-05-01
// Last changed: 2010-01-21


#include <dolfin.h>
#include <dolfin/fem/UFC.h>
#include <dolfin/main/MPI.h>

#include "unicorn/MeshBC.h"
#include "unicorn/NodeNormal.h"

#include <cstring>
#include <map>

#define max(a,b) (a > b ? a : b) ;

using namespace dolfin;
using namespace dolfin::unicorn;

//-----------------------------------------------------------------------------
MeshBC::MeshBC(Mesh& mesh, SubDomain& sub_domain, MeshFunction<bool>& cells,
	       GenericVector* node_values)
  : BoundaryCondition(), mesh(mesh),
    sub_domains(0), sub_domain(0), sub_domains_local(false),
    user_sub_domain(&sub_domain), cells(cells), As(0),
    row_block(0), zero_block(0), a1_indices_array(0),
    boundary(0), cell_map(0), vertex_map(0), node_values(node_values)
{
 // Initialize sub domain markers
  init(sub_domain);
  
  sub_system = SubSystem(0);
}
//-----------------------------------------------------------------------------
MeshBC::~MeshBC()
{
  // Delete sub domain markers if created locally
  if ( sub_domains_local )
    delete sub_domains;
  
  if( As )
    delete As;
  
  if(a1_indices_array)
    delete [] a1_indices_array; 
  if(row_block)
    delete [] row_block;
  if(zero_block)
    delete [] zero_block;
  if(boundary)
    delete boundary;
}
//-----------------------------------------------------------------------------
void MeshBC::apply(GenericMatrix& A, GenericVector& b, 
		   const Form& form)
{
  apply(A, b, form.dofMaps()[1], form);
}
//-----------------------------------------------------------------------------
void MeshBC::apply(GenericMatrix& A, GenericVector& b, const GenericVector& x, const Form& form)
{
  apply(A, b, form.dofMaps()[1], form);
}
//-----------------------------------------------------------------------------
void MeshBC::apply(GenericMatrix& A, GenericVector& b, const GenericVector& x, const DofMap& dof_map,
		   const ufc::form& form)
{
  apply(A, b, dof_map, form);
}
//-----------------------------------------------------------------------------
void MeshBC::apply(GenericMatrix& A, GenericVector& b, 
			const DofMap& dof_map, const ufc::form& form)
{

}
//-----------------------------------------------------------------------------
void MeshBC::apply(GenericMatrix& A, GenericVector& b, const DofMap& dof_map, 
		   const Form& form)
{
  
  message("Applying MeshBC boundary conditions to linear system.");
  //dolfin_set("output destination", "silent");
  
  UFC ufc(form.form(), mesh, form.dofMaps());

  if(boundary == 0)  {
    boundary = new BoundaryMesh(mesh);
    if (boundary->numCells()) 
    {
      cell_map = boundary->data().meshFunction("cell map");
      vertex_map = boundary->data().meshFunction("vertex map");
    }
  }    

  if(As == 0) {
        
    // Create data structure for local assembly data    

    const std::string la_backend = dolfin_get("linear algebra backend");
    if(la_backend == "JANPACK") {
      As = new Matrix(A.size(0), A.size(1));    
      *(As->instance()) =  A;
      //      (*(As->instance())).down_cast<JANPACKMat>().dup(A);
    }
    else {
      As = new Matrix();    
      (*(As->instance())).down_cast<PETScMatrix>().dup(A);
    }
    
    if(MPI::numProcesses() > 1) {
      std::map<uint, uint> mapping;
      for(CellIterator c(mesh); !c.end(); ++c) {
	ufc.update(*c, mesh.distdata());
	(form.dofMaps())[0].tabulate_dofs(ufc.dofs[0], ufc.cell, c->index());
	
	for(uint j = 0; j < (form.dofMaps())[0].local_dimension(); j++) 
	  off_proc_rows.insert(ufc.dofs[0][j]);
      }
      
      b.init_ghosted(off_proc_rows.size(), off_proc_rows, mapping);
    }

    row_block = new real[A.size(0)];
    zero_block = new real[A.size(0)]; 
    a1_indices_array = new uint[A.size(0)];
  }    
  
  // Copy global stiffness matrix into temporary one
  *(As->instance()) =  A;

  Array<uint> nodes;
  uint gdim = mesh.geometry().dim();
  uint cdim = mesh.type().numVertices(mesh.topology().dim());

  int d = mesh.topology().dim();

  uint count = 0;    
  for (CellIterator c(mesh); !c.end(); ++c) {
    if(cells(*c))
    {
    for (VertexIterator v(*c); !v.end(); ++v) {
      
      Vertex vertex = *v;
      
      uint node = vertex.index();     
      if(!mesh.distdata().is_ghost(node, 0) || MPI::numProcesses() == 1) {
	
	Cell cell(mesh, (vertex.entities(d))[0]);
	
	uint *cvi = cell.entities(0);
	uint ci = 0;
	for(ci = 0; ci < cell.numEntities(0); ci++)
	  if(cvi[ci] == node)
	    break;
	
	ufc.update(cell, mesh.distdata());
	(form.dofMaps())[0].tabulate_dofs(ufc.dofs[0], ufc.cell, cell.index());
	
	
	for(uint i = 0; i < gdim; i++, ci+=cdim) 
	  nodes.push_back(ufc.dofs[0][ci]); 
	
	
	applyMeshBC((Matrix&) A, *As, (Vector&) b, mesh, node, nodes); 
	count++;   
	nodes.clear();
      }
    }
    }
  }

  // Apply changes in the temporary matrix
  As->apply();


  // Apply changes in the stiffness matrix and load vector
  A = *(As->instance());
  b.apply();

}
//-----------------------------------------------------------------------------
void MeshBC::init(SubDomain& sub_domain)
{
  // Create mesh function for sub domain markers on facets
  mesh.init(0);
  sub_domains = new MeshFunction<uint>(mesh, 0);
  sub_domains_local = true;

  // Mark everything as sub domain 1
  (*sub_domains) = 1;
  
  // Mark the sub domain as sub domain 0
  sub_domain.mark(*sub_domains, 0);
}

//-----------------------------------------------------------------------------
void MeshBC::applyMeshBC(Matrix& A, Matrix& As, Vector& b, Mesh& mesh, 
			 uint node, Array<uint>& nodes)
{

  int nsdim = mesh.topology().dim();
  
  // Get number of nozero elements in the rows of A
  //FIXME implement nzmax in dolfin 0.8
  int nzm = A.size(0);

  // Now, we use node_type vector, which is defined for all nodes at the boundary:
  // node_type = 1: node at the surface; Apply slip BC as usual
  // node_type = 2: node at the edge; Apply slip BC in modified way using the secong normal which is tau_1
  // node_type = 3: node at the corner; Apply no slip BC, meening that u = 0 at this node
  
  // Note a2 is just index(a1) + N
  
  // define type of node, in 3D: 1 - surface; 2 - edge; 3 - corner; 
  //                      in 2D: 1 - surface; 2 - corner;
  
  uint a1_ncols, a2_ncols, a3_ncols; 
  Array<real> a1, a2, a3;
  Array<uint> a1_indices, a2_indices, a3_indices;
  a1_ncols = a2_ncols = a3_ncols = 0;

  A.getrow(nodes[0], a1_indices, a1); 
  a1_ncols = a1_indices.size();
  A.getrow(nodes[1], a2_indices, a2);

  a2_ncols = a2_indices.size();
  if (nsdim == 3) {
    A.getrow(nodes[2], a3_indices, a3);		  
    a3_ncols = a3_indices.size();
  }
  
  nzm = a1_ncols;
  memset(row_block, 0.0, nzm * sizeof(real));
  memset(zero_block, 0.0, nzm * sizeof(real));

  // Get RHS for the row node
  real l1, l2, l3;
  l3 = 0;
  l1 = b[nodes[0]];
  l2 = b[nodes[1]];
  if (nsdim == 3)
    l3 = b[nodes[2]];
    
  std::copy(a1_indices.begin(), a1_indices.end(), a1_indices_array);

  // Apply no-slip for the node on the corner  in 3D and 2D
  if (true) {
    As.set(zero_block, 1, &nodes[0], static_cast<uint>(a1_ncols), a1_indices_array);
    As.set(zero_block, 1, &nodes[1], static_cast<uint>(a1_ncols), a1_indices_array);
    if(nsdim == 3) 
      As.set(zero_block, 1, &nodes[2], static_cast<uint>(a1_ncols), a1_indices_array);
    
    Aset(As, nodes[0], nodes[0], 1.0);
    Aset(As, nodes[1], nodes[1], 1.0);
    if(nsdim == 3)
      Aset(As, nodes[2], nodes[2], 1.0);

    if(node_values)
    {
      bset(b, nodes[0], (*node_values)[nodes[0]]);
      bset(b, nodes[1], (*node_values)[nodes[1]]);
      if(nsdim == 3)
	bset(b, nodes[2], (*node_values)[nodes[2]]);
    }
    else
    {
      bset(b, nodes[0], 0.0);
      bset(b, nodes[1], 0.0);
      if(nsdim == 3)
	bset(b, nodes[2], 0.0);    
    }
  }
}
//-----------------------------------------------------------------------------

