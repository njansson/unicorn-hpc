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

#include "unicorn/SlipBC.h"
#include "unicorn/NodeNormal.h"
#include <map>

#define max(a,b) (a > b ? a : b) ;

using namespace dolfin;

//-----------------------------------------------------------------------------
SlipBC::SlipBC(Mesh& mesh, SubDomain& sub_domain, NodeNormal& Node_normal)
  : BoundaryCondition(), mesh(mesh),
    sub_domains(0), sub_domain(0), sub_domains_local(false),
    user_sub_domain(&sub_domain), node_normal(Node_normal), As(0),
    row_block(0), zero_block(0), a1_indices_array(0),
    boundary(0), cell_map(0), vertex_map(0)
{
 // Initialize sub domain markers
  init(sub_domain);
  
  sub_system = SubSystem(0);
}
//-----------------------------------------------------------------------------
SlipBC::SlipBC(MeshFunction<uint>& sub_domains,
	       uint sub_domain)
  : BoundaryCondition(), mesh(sub_domains.mesh()),
    sub_domains(&sub_domains), sub_domain(sub_domain), sub_domains_local(false),
    node_normal(mesh), As(0),  row_block(0), zero_block(0), a1_indices_array(0),
    boundary(0), cell_map(0), vertex_map(0)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
SlipBC::SlipBC(Mesh& mesh,
	       SubDomain& sub_domain,
	       const SubSystem& sub_system)
  : BoundaryCondition(), mesh(mesh),
    sub_domains(0), sub_domain(0), sub_domains_local(false),
    sub_system(sub_system), user_sub_domain(&sub_domain), node_normal(mesh),
    As(0),  row_block(0), zero_block(0), a1_indices_array(0),
    boundary(0), cell_map(0), vertex_map(0)
{
  // Set sub domain markers
  init(sub_domain);
  //Function g(0);
}

//-----------------------------------------------------------------------------
SlipBC::SlipBC(MeshFunction<uint>& sub_domains,
	       uint sub_domain,
	       const SubSystem& sub_system)
  : mesh(sub_domains.mesh()),
    sub_domains(&sub_domains), sub_domain(sub_domain), sub_domains_local(false),
    sub_system(sub_system), node_normal(node_normal), As(0), row_block(0),
    zero_block(0), a1_indices_array(0), boundary(0), cell_map(0), vertex_map(0)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
SlipBC::~SlipBC()
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
void SlipBC::apply(GenericMatrix& A, GenericVector& b, 
		   const Form& form)
{
  apply(A, b, form.dofMaps()[1], form);
}
//-----------------------------------------------------------------------------
void SlipBC::apply(GenericMatrix& A, GenericVector& b, const GenericVector& x, const Form& form)
{
  apply(A, b, form.dofMaps()[1], form);
}
//-----------------------------------------------------------------------------
void SlipBC::apply(GenericMatrix& A, GenericVector& b, const GenericVector& x, const DofMap& dof_map,
		   const ufc::form& form)
{
  apply(A, b, dof_map, form);
}
//-----------------------------------------------------------------------------
void SlipBC::apply(GenericMatrix& A, GenericVector& b, 
			const DofMap& dof_map, const ufc::form& form)
{

}
//-----------------------------------------------------------------------------
void SlipBC::apply(GenericMatrix& A, GenericVector& b, const DofMap& dof_map, 
		   const Form& form)
{
  
  if(MPI::processNumber() == 0)
    dolfin_set("output destination", "terminal");
  message("Applying SlipBC boundary conditions to linear system.");
  dolfin_set("output destination", "silent");
  
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
    As = new Matrix();
    (*(As->instance())).down_cast<PETScMatrix>().dup(A);
    
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

  uint count = 0;    
  if (boundary->numCells()) {
    for (VertexIterator v(*boundary); !v.end(); ++v) {
      
      Vertex vertex(mesh, vertex_map->get(*v));
      
      // Skip facets not inside the sub domain
      if ( (*sub_domains)(vertex) != sub_domain ) 
	continue;
      
      uint node = vertex.index();     
      if(!mesh.distdata().is_ghost(node, 0) || MPI::numProcesses() == 1) {
	
	Cell cell(mesh, (vertex.entities(3))[0]);
	
	uint *cvi = cell.entities(0);
	uint ci = 0;
	for(ci = 0; ci < cell.numEntities(0); ci++)
	  if(cvi[ci] == node)
	    break;
	
	ufc.update(cell, mesh.distdata());
	(form.dofMaps())[0].tabulate_dofs(ufc.dofs[0], ufc.cell, cell.index());
	
	
	for(uint i = 0; i < gdim; i++, ci+=cdim) 
	  nodes.push_back(ufc.dofs[0][ci]); 
	
	
	applySlipBC((Matrix&) A, *As, (Vector&) b, mesh, node, nodes); 
	count++;   
	nodes.clear();
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
void SlipBC::init(SubDomain& sub_domain)
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
void SlipBC::applySlipBC(Matrix& A, Matrix& As, Vector& b, Mesh& mesh, 
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
  int n_type = node_normal.node_type.get(node);


  // Apply no-slip for the node on the corner  in 3D and 2D
  if ((n_type == 3 && nsdim == 3) || (n_type == 2 && nsdim == 2)) {
    As.set(zero_block, 1, &nodes[0], static_cast<uint>(a1_ncols), a1_indices_array);
    As.set(zero_block, 1, &nodes[1], static_cast<uint>(a1_ncols), a1_indices_array);
    if(nsdim == 3) 
      As.set(zero_block, 1, &nodes[2], static_cast<uint>(a1_ncols), a1_indices_array);
    
    Aset(As, nodes[0], nodes[0], 1.0);
    Aset(As, nodes[1], nodes[1], 1.0);
    if(nsdim == 3)
      Aset(As, nodes[2], nodes[2], 1.0);

    bset(b, nodes[0], 0.0);
    bset(b, nodes[1], 0.0);
    if(nsdim == 3)
      bset(b, nodes[2], 0.0);    
  }
  else {

    // get components of tau
    real t11, t12, t13;
    real t21, t22, t23;
    t11 = t12 = t21 = t22 = t13 = t23 = 0;
    if (nsdim == 3) {
      t11 = node_normal.tau_1[0].get(node);
      t12 = node_normal.tau_1[1].get(node);
      t13 = node_normal.tau_1[2].get(node);

      t21 = node_normal.tau_2[0].get(node);
      t22 = node_normal.tau_2[1].get(node);
      t23 = node_normal.tau_2[2].get(node);
    }

    // get componets of normal
    real n1, n2, n3;
    n3 = 0;
    n1 = node_normal.normal[0].get(node);
    n2 = node_normal.normal[1].get(node);

    if (nsdim == 3)
      n3 = node_normal.normal[2].get(node);

    // find maximum of the normal components and put them to the diagonal 
    real maxn = max(fabs(n1), fabs(n2));
    if (nsdim == 3)
      maxn = max(maxn, fabs(n3));

    // r1, r2, r3 are rows which are corresponds to the boundary point
    // find the maximum component of the normal and put it to the diagonal
    uint r1 = 0; 
    uint r2 = 0; 
    uint r3 = nodes[2]; 
      
    if (fabs(fabs(n1) - maxn) < DOLFIN_EPS) {       // n1 is a largest component
      r1 = nodes[0];
      r2 = nodes[1];
    }
    else if (fabs(fabs(n2) - maxn) < DOLFIN_EPS) {  // n2 is a largest component
      r1 = nodes[1];
      r2 = nodes[0];
    }
    else {                         // n3 is a largest component
      if (nsdim == 3){
	r1 = nodes[2];
	r2 = nodes[1]; r3 = nodes[0];
      }
    }
      
    As.set(zero_block, 1, &r1, static_cast<uint>(a1_ncols), a1_indices_array);

    Aset(As, r1, nodes[0], n1);
    Aset(As, r1, nodes[1], n2);      
    if (nsdim == 3)
      Aset(As, r1, nodes[2], n3);

    bset(b, r1, 0.0);

    // set the row using tangent vector in 2D
    /*
    if (nsdim == 2) {
      rows[0] = r2;
      
      // Set rows of matrix to the identity matrix
      for(uint i = 0; i < a1_ncols; i++)
	row_block[i] = (a1[i] * t1 + a2[i] * t2);
      
      As.set(row_block, 1, rows, static_cast<uint>(a1_ncols), a1_indices_array);
      
      bset(b, r2, (l1 * t1 + l2 * t2));
    }
    */

    // if node is on the edge in 3D case:
    if (n_type == 2 && nsdim == 3) {
      
      // The variables below are used in rearanging rows
      real nn1, nn2, nn3, tt1, tt2, tt3;
      
      int non0 = 0;
      if (fabs(n1)<DOLFIN_EPS) non0++;
      if (fabs(n2)<DOLFIN_EPS) non0++;
      if (fabs(n3)<DOLFIN_EPS) non0++;
      
      if (non0 > 1) {
	nn1 = n1; nn2 = n2; nn3 = n3;
	tt1 = t11; tt2 = t12; tt3 = t13;
      }
      else {
	tt1 = n1; tt2 = n2; tt3 = n3;
	nn1 = t11; nn2 = t12; nn3 = t13;
      }
      
      // find maximum of the normal components and put them to the diagonal 
      real maxn = max(fabs(nn1), fabs(nn2));
      maxn = max(maxn, fabs(nn3));
      
      // r1, r2, r3 are rows which are corresponds to the boundary point
      // find the maximum component of the normal and put it to the diagonal
      r1 = 0; 
      r2 = 0; 
      r3 = nodes[2];
      
      if (fabs(fabs(nn1) - maxn) < DOLFIN_EPS) {       // n1 is a largest component
	r1 = nodes[0];
	r2 = nodes[1];
      }
      else if (fabs(fabs(nn2) - maxn) < DOLFIN_EPS) {   // n2 is a largest component
	r1 = nodes[1];
	r2 = nodes[0];
      }
      else {                         // n3 is a largest component
	r1 = nodes[2];
	r2 = nodes[1]; 
	r3 = nodes[0];
      }
      

      As.set(zero_block, 1, &r1, static_cast<uint>(a1_ncols), a1_indices_array);

      Aset(As, r1, nodes[0], nn1);
      Aset(As, r1, nodes[1], nn2);
      Aset(As, r1, nodes[2], nn3);
	  

      maxn = 0.0;
      maxn = max(fabs(tt1), fabs(tt2)); 
      maxn = max(maxn, fabs(tt3));
	  
      // Start rearanging rows
      if (fabs(fabs(tt1) - maxn) < DOLFIN_EPS) {       // tt1 is a largest component
	r2 = nodes[0]; 
	r3 = (r1 == nodes[1] ?  nodes[2] : nodes[1]) ;
	if (r2 == r1) {
	  real mm = max(fabs(tt2), fabs(tt3));
	  if ( fabs(fabs(tt2) - mm) < DOLFIN_EPS) {
	    r2 = nodes[1]; 
	    r3 = nodes[2];
	  }
	  else {
	    r2 = nodes[2]; 
	    r3 = nodes[1];
	  }
	}
      }
      else if (fabs(fabs(tt2) - maxn) < DOLFIN_EPS) {  // tt2 is a largest component
	r2 = nodes[1];
	r3 = (r1 == nodes[0] ? nodes[2] : nodes[0]) ;
	if (r2 == r1) {
	  real mm = max(fabs(tt1), fabs(tt3));
	  if (fabs(fabs(tt1) - mm) < DOLFIN_EPS) {
	    r2 = nodes[0]; 
	    r3 = nodes[2];
	  }
	  else {
	    r2 = nodes[2]; 
	    r3 = nodes[0];
	  }
	}
      }
      else {                        // tt3 is a largest component
	r2 = nodes[2]; 
	r3 = (r1 == nodes[0] ? nodes[1] : nodes[0]) ;
	if (r2 == r1) {
	  real mm = max(fabs(tt1), fabs(tt2));
	  if (fabs(fabs(tt1) - mm) < DOLFIN_EPS) {
	    r2 = nodes[0]; 
	    r3 = nodes[1];
	  }
	  else {
	    r2 = nodes[1]; 
	    r3 = nodes[0];
	  }
	}
      }
      
      As.set(zero_block, 1, &r2, static_cast<uint>(a1_ncols), a1_indices_array);

      Aset(As, r2, nodes[0], tt1);
      Aset(As, r2, nodes[1], tt2);
      Aset(As, r2, nodes[2], tt3);
      bset(b, r2, 0.0);

	  
      // Set second row of slip BC (tau..) to the row r3
	  
      // Set rows of matrix to the identity matrix
      for(uint i = 0; i < a1_ncols; i++)
	row_block[i] = a1[i] * t21 + a2[i] * t22 + a3[i] * t23;
	  
      //As.set(row_block, rows, 1, a1_indices_array, a1_ncols);
      As.set(row_block, 1, &r3, static_cast<uint>(a1_ncols), a1_indices_array);
      bset(b, r3, l1 * t21 + l2 * t22 + l3 * t23);
	  
    }

    // The case where node is on the surface
    if (n_type == 1 && nsdim == 3) {	      
      int ind1 = 0;
      int ind2 = 0;
      
      for (uint i = 0; i < a1_ncols; i++) {
	if (r2 == a1_indices[i])
	  ind1 = i;
	if (r3 == a1_indices[i])
	  ind2 = i;
      }

      // find maximum of the elements of 2 vectors which are belong to diagonal
      real row_d11 = a1[ind1] * t11 + a2[ind1] * t12 + a3[ind1] * t13 ;
      real row_d21 = a1[ind1] * t21 + a2[ind1] * t22 + a3[ind1] * t23 ;
      real row_d12 = a1[ind2] * t11 + a2[ind2] * t12 + a3[ind2] * t13 ;
      real row_d22 = a1[ind2] * t21 + a2[ind2] * t22 + a3[ind2] * t23 ;

      real maxr = max(fabs(row_d11), fabs(row_d21));
      maxr = max(maxr, fabs(row_d12));
      maxr = max(maxr, fabs(row_d22));
	  
      // define new rows according the above maximum:
      uint row_d1 = 0;
      uint row_d2 = 0;
      if (fabs( fabs(row_d11) - maxr) < DOLFIN_EPS  || fabs( fabs(row_d22) - maxr) < DOLFIN_EPS) {
	row_d1 = r2;
	row_d2 = r3;
      }
      else {
	row_d1 = r3;
	row_d2 = r2;
      }
	  
      for(uint i = 0; i < a1_ncols; i++)
	row_block[i] = a1[i] * t11 + a2[i] * t12 + a3[i] * t13;

      //As.set(row_block, rows, 1, a1_indices_array, a1_ncols);
      As.set(row_block, 1, &row_d1, static_cast<uint>(a1_ncols), a1_indices_array);
      bset(b, row_d1, l1 * t11 + l2 * t12 + l3 * t13);

      // Set second row of slip BC (tau..) to the row r3
      for(uint i = 0; i < a1_ncols; i++)
	row_block[i] = a1[i] * t21 + a2[i] * t22 + a3[i] * t23;
	  
      //As.set(row_block, rows, 1, a1_indices_array, a1_ncols);
      As.set(row_block, 1, &row_d2, static_cast<uint>(a1_ncols), a1_indices_array);
      bset(b, row_d2, l1 * t21 + l2 * t22 + l3 * t23);

    }    
  }
}
//-----------------------------------------------------------------------------

