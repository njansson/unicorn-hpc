// Copyright (C) 2007 Murtazo Nazarov
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2008-2009.
//
// First added:  2007-05-01
// Last changed: 2009-12-30
                                                                                         
#include <dolfin.h>
#include <dolfin/main/MPI.h>
#include "unicorn/NodeNormal.h"

#include <map>

#define B(row,col,nrow) ((row) + ((nrow)*(col)))

using namespace dolfin;
using namespace dolfin::unicorn;

//-----------------------------------------------------------------------------
NodeNormal::NodeNormal(NodeNormal& node_normal) : mesh(node_normal.mesh),
						  normal(0), 
						  tau(0), tau_1(0), tau_2(0)
{
  *this = node_normal;
}
//-----------------------------------------------------------------------------
NodeNormal::NodeNormal(Mesh& mesh):mesh(mesh), normal(0), 
				   tau(0), tau_1(0), tau_2(0)
{

  normal = new MeshFunction<real>[mesh.topology().dim()];
  tau = new MeshFunction<real>[mesh.topology().dim()];
  tau_1 = new MeshFunction<real>[mesh.topology().dim()];
  tau_2 = new MeshFunction<real>[mesh.topology().dim()];

  for(uint i = 0; i < mesh.topology().dim(); i++) {
    normal[i].init(mesh, 0);
    tau[i].init(mesh, 0);
    tau_1[i].init(mesh, 0);
    tau_2[i].init(mesh, 0);
  }
  node_type.init(mesh, 0);
 
  ComputeNormal(mesh);
}
//-----------------------------------------------------------------------------
NodeNormal::~NodeNormal()
{
  clear();
}
//-----------------------------------------------------------------------------
void NodeNormal::clear()
{
  if ( normal )
    delete[] normal;
  if ( tau )
    delete[] tau;
  if ( tau_1 )
    delete[] tau_1;
  if ( tau_2 )
    delete[] tau_2;  
}
//-----------------------------------------------------------------------------
NodeNormal& NodeNormal::operator=(NodeNormal& node_normal)
{

  clear();
  
  normal = new MeshFunction<real>[mesh.topology().dim()];
  tau = new MeshFunction<real>[mesh.topology().dim()];
  tau_1 = new MeshFunction<real>[mesh.topology().dim()];
  tau_2 = new MeshFunction<real>[mesh.topology().dim()];
  
  for(uint i = 0; i < mesh.topology().dim(); i++) {
    normal[i].init(mesh, 0);
    tau[i].init(mesh, 0);
    tau_1[i].init(mesh, 0);
    tau_2[i].init(mesh, 0);

    for(VertexIterator v(mesh); !v.end(); ++v)  
    {
      normal[i].set(*v, node_normal.normal[i].get(*v));
      tau[i].set(*v, node_normal.tau[i].get(*v));
      tau_1[i].set(*v, node_normal.tau_1[i].get(*v));
      tau_2[i].set(*v, node_normal.tau_2[i].get(*v));    
    }
      
  }
  node_type.init(mesh, 0);
  for(VertexIterator v(mesh); !v.end(); ++v)
    node_type.set(*v, node_normal.node_type.get(*v));

  return *this;
}
//-----------------------------------------------------------------------------
void NodeNormal::ComputeNormal(Mesh& mesh)
{

  mesh.renumber();
  uint rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();
  Array<real> *send_buff_type = new Array<real>[pe_size];
  Array<uint> *send_buff_index = new Array<uint>[pe_size];
  _map<uint, bool> used_shared;

  for(MeshSharedIterator s(mesh.distdata(), 0); !s.end(); ++s)
    used_shared[s.index()] = false;

  int nsdim = mesh.topology().dim();
  
  // Iterate over all cells in the boundary mesh
  BoundaryMesh boundary(mesh);
  MeshFunction<uint>* cell_map = boundary.data().meshFunction("cell map");
  MeshFunction<uint>* vertex_map = boundary.data().meshFunction("vertex map");
  if(boundary.numCells() > 0) 
    Progress p_boundary("Assembling boundary contributions", boundary.numCells());
  
  real normal_vec[3];
  real normal_vec1[3];
  real tau_vec[2];
  real tau_vec_1[3];
  real tau_vec_2[3];

  real *n_block = 0;
  real *ns1_block = 0;
  real *ns2_block = 0;
  real *area_block = 0;
  
  
  if(MPI::numProcesses() > 1)
    cache_shared_area(mesh, boundary,nsdim, vertex_map, cell_map);
  
  // Computation of normals to the boundary vertices
  if(boundary.numCells() > 0) {
    for (VertexIterator n(boundary); !n.end(); ++n) {
      
      int id = vertex_map->get(*n);
      
      int Nrow = 0;
      for (CellIterator boundary_cell(*n); !boundary_cell.end(); ++boundary_cell)
	Nrow++;
      
      // Add storage for shared vertices cell normals      
      if(MPI::numProcesses() > 1) 
	if(mesh.distdata().is_shared(id, 0))
	  Nrow = num_cells[mesh.distdata().get_global(id, 0)];
      
      int Ncol = 3;
      
      n_block = new real[Nrow * Ncol]; // contains of all normals from the surfaces 
      
      int irow = 0;
      area_block = new real[Nrow]; // contains of all areas from the elements 
      
      real sum_area = 0.0;         // sum of all areas of the elements
      
      // in a case of node on the edge
      real sum_area_s1 = 0.0;      // sum of all areas of the elements which belong to S1
      real sum_area_s2 = 0.0;      // sum of all areas of the elements which belong to S2
      
      for (int i = 0; i < nsdim; i++) {
	normal_vec[i] = 0.0;
	normal_vec1[i] = 0.0;
      }
      
      //for (CellIterator boundary_cell(boundary); !boundary_cell.end(); ++boundary_cell)
      for (CellIterator boundary_cell(*n); !boundary_cell.end(); ++boundary_cell) 
      {
	// Create mesh facet corresponding to boundary cell
	Facet mesh_facet(mesh, cell_map->get(*boundary_cell));
	dolfin_assert(mesh_facet.numEntities(nsdim) == 1);
	
	// Get cell to which facet belongs (pick first, there is only one)
	Cell mesh_cell(mesh, mesh_facet.entities(mesh.topology().dim())[0]);	  
	// Get local index of facet with respect to the cell                      
	uint local_facet = mesh_cell.index(mesh_facet);
	
	// Compute the lenght/area of the cell
	real area = mesh_cell.volume();
	area_block[irow] = area;
	// Get sum of the length/area of all neighbouring cells
	sum_area += area;
	
	// Compute a normal to the boundary
	normal_vec[0] += area * mesh_cell.normal(local_facet, 0);
	normal_vec[1] += area * mesh_cell.normal(local_facet, 1);
	if (nsdim == 3)
	  normal_vec[2] += area * mesh_cell.normal(local_facet, 2);
	
	for (int j = 0; j < Ncol; j++)
	  n_block[B(irow, j, Nrow)] = mesh_cell.normal(local_facet, j);
	irow ++;
      }
      

      if(MPI::numProcesses() > 1)
	if(mesh.distdata().is_shared(id, 0)) {
	  uint glb_index = mesh.distdata().get_global(id, 0);
	  //	  sum_area += shared_area[glb_index];
	  normal_vec[0] += shared_normal[shared_noffset[glb_index]];
	  normal_vec[1] += shared_normal[shared_noffset[glb_index]+1];
	  if(nsdim == 3)
	    normal_vec[2] += shared_normal[shared_noffset[glb_index]+2];
	  
	  
	  Array<real> tmp = shared_area_block[glb_index];
	  std::vector<real>::iterator iter;
	  
	  int iirow = irow;
	  for(iter = tmp.begin(); iter != tmp.end(); iter++){
	    area_block[iirow] = *iter;
	    sum_area += area_block[iirow];
	    iirow++;
	  }
	  
	  tmp = normal_block[glb_index];
	  for(iter = tmp.begin(); iter != tmp.end(); iter +=3){
	    for(int j = 0; j < Ncol; j++)
	      n_block[B(irow,j,Nrow)] = *(iter + j);
	    irow++;	    
	  }
	}
      

      ns1_block = new real[Nrow * Ncol]; // contains of all normals from the case S1
      ns2_block = new real[Nrow * Ncol]; // contains of all normals from the case S2

      real n1[3];
      int n_type = 1;                 // n_type is a node or vertex type: 
                                      // n_type = 1 vector in S1, 
                                      // n_type = 2 in S2, 
                                      // n_type = 3 in S3
                                      // S1 - vertex lies on the surface, 
                                      // S2 - edge, S3 - corner
      real alpha = 0.0;               // needed for the angle between the vectors

      // take one normal from the matrix of normals n1_block
      for (int j = 0; j < Ncol; j++)
	n1[j] =  n_block[B(0,j,Nrow)];

      if(!ParameterSystem::parameters.defined("alpha_max"))
	dolfin_add("alpha_max", DOLFIN_PI/3.0);
      real alpha_max = (real) dolfin_get("alpha_max");
      int Nrow_ns1 = 0; // it will be number of rows for the ns1_block
      int Nrow_ns2 = 0; // it will be number of rows for the ns2_block

      // find all vector which are from S2 - edge
      
      for (int i = 0; i < Nrow; i++) {
	real csum1 = 0.0;
	real csum2 = 0.0;
	
	for (int in = 0; in < nsdim; in ++) {
	  csum1 += sqr(n1[in]);
	  csum2 += sqr(n_block[B(i, in, Nrow)]);
	}

	real n_norm  = sqrt(csum1);
	real nj_norm = sqrt(csum2);
	
	real csum3 = 0.0;
	for (int in = 0; in < nsdim; in ++)
	  csum3 += n1[in] * n_block[B(i,in,Nrow)];
	
	alpha = acos ( csum3 / (n_norm * nj_norm ));

	if (nsdim == 2 && fabs(alpha) > alpha_max)
	  n_type = 2;
	  
	if (nsdim == 3 ) {
	  
	  if (fabs(alpha) > alpha_max) {
	    n_type = 2;
	    
	    // save normals to the block ns2_block
	    for (int j = 0; j < Ncol; j++) {
	      ns2_block[B(Nrow_ns2,j,Nrow)] = n_block[B(i,j,Nrow)];
	    }
	    
	    // calculate the area of all elements belongs ot S2
	    sum_area_s2 += area_block[i]; 
	    Nrow_ns2++;
	  }
	  else {
	    // save normals to the block ns1_block
	    for (int j = 0; j < Ncol; j++)
	      ns1_block[B(Nrow_ns1, j,Nrow)] =  n_block[B(i,j,Nrow)];
	    
	    // calculate the area of all elements belongs ot S2
	    sum_area_s1 += area_block[i]; 
	    Nrow_ns1++;	    
	  }
	}
      } // end of for

      if (nsdim == 3) {
	// go through all normals from S2 and find normals from S3 if any
	if (n_type == 2) {
	  for (int j = 0; j < Ncol; j++)
	    n1[j] =  ns2_block[B(0,j,Nrow)];
	  
	  // find all vector which are from S3 - corner
	  for (int i = 0; i < Nrow_ns2; i++) {
	    real n_norm = sqrt( sqr(n1[0]) + sqr(n1[1]) + sqr(n1[2]) );
	    real nj_norm = sqrt( sqr(ns2_block[B(i,0,Nrow)]) + 
				 sqr(ns2_block[B(i,1,Nrow)]) +
				 sqr(ns2_block[B(i,2,Nrow)]) );
	    alpha = acos ( (n1[0] * ns2_block[B(i,0,Nrow)] + n1[1] * 
			    ns2_block[B(i,1,Nrow)] + n1[2] * 
			    ns2_block[B(i,2,Nrow)] ) / n_norm * nj_norm );
	    if (fabs(alpha) > alpha_max) {
	      n_type = 3;
	    }
	  }
	}
	
	// Compute a normal to the vertex
	if (n_type == 1) { // the vertex lies on the surface, case S1
	  for (int in = 0; in < nsdim; in++) {
	    normal_vec[in] /= sum_area; 
	  }
	}
	else if (n_type == 2) { // the vertex lies on the edge, case S2
	  // now we use the matrix of normals for S2: ns1_block and ns2_block
	  // note ns1_block is used for the first normal and 
	  // ns2_block is for the second normal at the same edge
	  // here we take the second normal instead of one tangental vector
	  
	  // compute the first normal from the normals of S1
	  for (int in = 0; in < nsdim; in++)
	    normal_vec[in] = 0.0;
	  
	  for (int i = 0; i < Nrow_ns1; i++)
	    for (int in = 0; in < nsdim; in++)
	      normal_vec[in] += ns1_block[B(i,in,Nrow)];
	  
	  for (int in = 0; in < nsdim; in++)
	    normal_vec[in] /= sum_area_s1; 
	  
	  // compute the second normal from the normals of S2
	  for (int i = 0; i < Nrow_ns2; i++)
	    for (int in = 0; in < nsdim; in++)
	      normal_vec1[in] += ns2_block[B(i,in,Nrow)];
	  
	  for (int in = 0; in < nsdim; in++)
	    normal_vec1[in] /= sum_area_s2; 
	}
      }
      
      // Find the length of the normal vector and normal_vecize it
      real l = 0.0; 
      for (int in = 0; in < nsdim; in++)
	l += sqr(normal_vec[in]);
      l = sqrt(l);
      
      for (int in = 0; in < nsdim; in++) {
	normal_vec[in] /= l;
	
	if (fabs(normal_vec[in]) < DOLFIN_EPS) // check if the component is too small replace it by zero
	  normal_vec[in] = 0.0;
      }
      
      if (nsdim == 3 && n_type == 2) { // normalize normal1 if the case S2 exist
	l = 0.0; 
	for (int in = 0; in < nsdim; in++)
	  l += sqr(normal_vec1[in]);
	l = sqrt(l);
	
	for (int in = 0; in < nsdim; in++) {
	  normal_vec1[in] /= l;
	  
	  if (fabs(normal_vec1[in]) < DOLFIN_EPS) // check if the component is too small replace it by zero
	    normal_vec1[in] = 0.0;
	}
	
      }
      
      if(!mesh.distdata().is_ghost(id, 0)) {      
	for (int in = 0; in < nsdim; in++)
	  normal[in].set(id, normal_vec[in]);
	node_type.set(id, n_type);
	used_shared[id] = true;
      }
      else{
	uint owner =  mesh.distdata().get_owner(id, 0);
	send_buff_type[ owner ].push_back( n_type );
	for (int in = 0; in < nsdim; in++)
	  send_buff_type[owner].push_back(normal_vec[in]);
	send_buff_index[ owner ].push_back(mesh.distdata().get_global(id, 0));
      }
      
      if (nsdim == 2) {
	error("nej");
	// 2D case
	tau_vec[0] =  normal_vec[1];
	tau_vec[1] = -normal_vec[0];
	
	// Find the length of the tau_vec vector and normalize it
	l = 0.0; 
	for (int in = 0; in < nsdim; in++)
	  l += sqr(tau_vec[in]);
	l = sqrt(l);
	
	for (int in = 0; in < nsdim; in++) {
	  tau_vec[in] /= l;
	
	  if (fabs(tau_vec[in]) < DOLFIN_EPS) // check if the component is too small replace it by zero
	    tau_vec[in] = 0.0;
	}	  
	if(!mesh.distdata().is_ghost(id, 0)){
	  error("nej nej");
	  /*
	  tau_vec_arr[0 * N + id] =  tau_vec[0];	 
	  tau_vec_arr[1 * N + id] = tau_vec[1];
	  */
	  //	  tau_vec_arr[0 * N + l2vidx[id]] =  tau_vec[0];	 
	  //	  tau_vec_arr[1 * N + l2vidx[id]] = tau_vec[1];
	}
	else {
	  uint owner =  mesh.distdata().get_owner(id, 0);
	  send_buff_type[owner].push_back(tau_vec[0]);
	  send_buff_type[owner].push_back(tau_vec[1]);
	}
      }
      else  {
	// 3D case
	real norm_inv = 0.0;
	if (fabs(normal_vec[0]) >= 0.5 || fabs(normal_vec[1]) >= 0.5) {
	  norm_inv = 1/sqrt(sqr(normal_vec[0] + sqr(normal_vec[1])));
	  if (n_type == 2) {
	    tau_vec_1[0] =  normal_vec1[0];
	    tau_vec_1[1] =  normal_vec1[1];
	    tau_vec_1[2] =  normal_vec1[2];	
	    
	  // cross product of tau_1 and normal_vec
	    tau_vec_2[0] = normal_vec[1] * tau_vec_1[2] - tau_vec_1[1] * normal_vec[2];
	    tau_vec_2[1] = normal_vec[2] * tau_vec_1[0] - tau_vec_1[2] * normal_vec[0];
	    tau_vec_2[2] = normal_vec[0] * tau_vec_1[1] - tau_vec_1[0] * normal_vec[1];
	  } 
	  else {
	    tau_vec_1[0] =  normal_vec[1] * norm_inv;
	    tau_vec_1[1] = -normal_vec[0] * norm_inv;
	    tau_vec_1[2] =  0.0;
	    
	    tau_vec_2[0] = -tau_vec_1[1] * normal_vec[2];
	    tau_vec_2[1] =  tau_vec_1[0] * normal_vec[2];
	    tau_vec_2[2] =  tau_vec_1[1] * normal_vec[0] - tau_vec_1[0] * normal_vec[1];
	  }
	  norm_inv = 0.0;
	}
	else { 
	  norm_inv = 1/sqrt(sqr(normal_vec[1] + sqr(normal_vec[2])));
	  if (n_type == 2) {
	    tau_vec_1[0] =  normal_vec1[0];
	    tau_vec_1[1] =  normal_vec1[1];
	    tau_vec_1[2] =  normal_vec1[2];
	  
	    // cross product of tau_vec_1 and normal_vec
	    tau_vec_2[0] = normal_vec[1] * tau_vec_1[2] - tau_vec_1[1] * normal_vec[2];
	    tau_vec_2[1] = normal_vec[2] * tau_vec_1[0] - tau_vec_1[2] * normal_vec[0];
	    tau_vec_2[2] = normal_vec[0] * tau_vec_1[1] - tau_vec_1[0] * normal_vec[1];
	  }
	  else {
	    tau_vec_1[0] =  0.0;
	    tau_vec_1[1] = -normal_vec[2] * norm_inv;
	    tau_vec_1[2] =  normal_vec[1] * norm_inv;	
	    
	    tau_vec_2[0] =  tau_vec_1[2] * normal_vec[1] - tau_vec_1[1] * normal_vec[2];
	    tau_vec_2[1] = -tau_vec_1[2] * normal_vec[0];
	    tau_vec_2[2] =  tau_vec_1[1] * normal_vec[0];
	  }
	  norm_inv = 0.0;
	}
	
	// Find the length of the tau_vec_1 vector and normalize it
	l = 0.0; 
	for (int in = 0; in < nsdim; in++)
	  l += sqr(tau_vec_1[in]);
	l = sqrt(l);
	for (int in = 0; in < nsdim; in++) {
	  tau_vec_1[in] /= l;
	  
	  // check if the component is too small replace it by zero
	  if (fabs(tau_vec_1[in]) < DOLFIN_EPS) 
	    tau_vec_1[in] = 0.0;
	}
	
	// Find the length of the tau_vec_2 vector and normalize it
	l = 0.0; 
	for (int in = 0; in < nsdim; in++) {
	  l += sqr(tau_vec_2[in]);
	}
	l = sqrt(l);
	for (int in = 0; in < nsdim; in++) {
	  tau_vec_2[in] /= l;
	  
	  if (fabs(tau_vec_2[in]) < DOLFIN_EPS) // check if the component is too small replace it by zero
	    tau_vec_2[in] = 0.0;
	}
	
	if(!mesh.distdata().is_ghost(id, 0)){

	  for (int in = 0; in < nsdim; in++) {
	    tau_1[in].set(id, tau_vec_1[in]);
	    tau_2[in].set(id, tau_vec_2[in]);
	  }	  
	}
	else {
	  uint owner =  mesh.distdata().get_owner(id, 0);
	  send_buff_type[owner].push_back(tau_vec_1[0]);
	  send_buff_type[owner].push_back(tau_vec_1[1]);
	  send_buff_type[owner].push_back(tau_vec_1[2]);
	  
	  send_buff_type[owner].push_back(tau_vec_2[0]);
	  send_buff_type[owner].push_back(tau_vec_2[1]);
	  send_buff_type[owner].push_back(tau_vec_2[2]);
	}
	
      }
      
      delete[] area_block;
      delete[] n_block;
      delete[] ns1_block;
      delete[] ns2_block;

    }    
  }

  if(MPI::numProcesses() > 1) {
    MPI_Status status;
    uint src,dest;
    int recv_count, recv_size, send_size, recv_count_data;
    recv_count_data = 0; // == 0
    for(uint i=0; i<pe_size; i++){
      send_size = send_buff_index[i].size();
      MPI_Reduce(&send_size, &recv_count, 1, MPI_INT,
		 MPI_SUM, i, dolfin::MPI::DOLFIN_COMM);
      
      send_size = send_buff_type[i].size();
      MPI_Reduce(&send_size, &recv_count_data, 1, MPI_INT,
		 MPI_SUM, i, dolfin::MPI::DOLFIN_COMM);
    }
    
    uint n_tau = nsdim + nsdim;
    if(nsdim == 3)
      n_tau = nsdim * 2 + nsdim;
    uint *recv_index =  new uint[recv_count];
    real *recv_type = new real[recv_count_data];
    
    
    for(uint i=1; i<pe_size; i++){
      src = (rank - i + pe_size) % pe_size;
      dest = (rank + i) % pe_size;
      
      MPI_Sendrecv(&send_buff_index[dest][0], send_buff_index[dest].size(),
		   MPI_UNSIGNED, dest, 0, recv_index, recv_count, MPI_UNSIGNED,
		   src, 0, dolfin::MPI::DOLFIN_COMM, &status);
      
      MPI_Sendrecv(&send_buff_type[dest][0], send_buff_type[dest].size(),
		   MPI_DOUBLE, dest, 1, recv_type, recv_count_data,
		   MPI_DOUBLE, src, 1, dolfin::MPI::DOLFIN_COMM, &status);
      MPI_Get_count(&status, MPI_DOUBLE, &recv_size);
      // Insert check if value assigned
      uint idx = 0;
      for(int j = 0; j<recv_size; j += (n_tau +1)){
	uint index = mesh.distdata().get_local(recv_index[idx], 0);
	if(!used_shared[index]) {
	  node_type.set(index, (int) recv_type[j]);
	  //	  node_vec_arr[l2vidx[index]] =  recv_type[j];
	  uint offset = 0;
	  if(nsdim == 2) {
	    error("nej");
	    for (int in = 0; in < nsdim; in++){
	      offset = j + in + 1;
	      //  normal_vec_arr[in * N + l2vidx[index]] = recv_type[offset];
	      //	      tau_vec_arr[in * N + l2vidx[index]] =  recv_type[offset+nsdim];
	    }
	  }
	  else {
	    for (int in = 0; in < nsdim; in++){
	      offset = j + in + 1;
	      normal[in].set(index, recv_type[offset]);
	      tau_1[in].set(index, recv_type[offset+nsdim]);
	      tau_2[in].set(index, recv_type[offset+2*nsdim]);
	    }
	  }
	  used_shared[index] = true;
	}
	idx++;
      }
    }

    delete[] recv_index;
    delete[] recv_type;
  }

  for(uint i=0; i<pe_size; i++){
    send_buff_type[i].clear();
    send_buff_index[i].clear();
  }
  delete[] send_buff_type;
  delete[] send_buff_index;


  shared_noffset.clear();
  num_cells.clear();
  shared_normal.clear();
  
  for (std::map<uint, Array<real> >::iterator it = normal_block.begin();
       it != normal_block.end(); it++) 
    it->second.clear();
  for (std::map<uint, Array<real> >::iterator it = shared_area_block.begin();
       it != shared_area_block.end(); it++) 
    it->second.clear();
  

  normal_block.clear();
  shared_area_block.clear();
}
//-----------------------------------------------------------------------------
void NodeNormal::cache_shared_area(Mesh& mesh, BoundaryMesh& boundary, uint nsdim, MeshFunction<uint> *vertex_map, MeshFunction<uint> *cell_map) 
{
  int rank = MPI::processNumber();
  int pe_size = MPI::numProcesses();

  // Send buff for global indices
  Array<uint> send_buff_indices;
  // Send buff for area and normals contribution,
  // stored as area, normal(0..2),
  Array<real> send_buff;
  Array<real> send_buffn;
  Array<uint> send_ni;
  Array<real> send_area;

  std::map<uint,bool> glb_boundary;

  _offset = 0;
  uint _noffset = 0;
  uint area_offset = 0;
  real normal[3];

  // Computation of normals to the boundary vertices
  if(boundary.numCells() > 0) {
    for (VertexIterator n(boundary); !n.end(); ++n) {   
      
      int id = vertex_map->get(*n);      
      
      glb_boundary[mesh.distdata().get_global(id, 0)] = true;
      
      // Cache number of cells for each shared vertex
      int Nrow = 0;
      if(mesh.distdata().is_shared(id, 0)) {	
	for (CellIterator boundary_cell(*n); 
	     !boundary_cell.end(); 
	     ++boundary_cell)  {
	  Nrow++;
	}
	num_cells[mesh.distdata().get_global(id, 0)] = Nrow;
      }
      int Ncol = 3;
      for (uint i = 0; i < nsdim; i++)
	normal[i] = 0.0;
      
      for (CellIterator boundary_cell(*n); !boundary_cell.end(); ++boundary_cell)
	{
	  // Create mesh facet corresponding to boundary cell
	  Facet mesh_facet(mesh, cell_map->get(*boundary_cell));
	  dolfin_assert(mesh_facet.numEntities(nsdim) == 1);
	  
	  // Get cell to which facet belongs (pick first, there is only one)
	  Cell mesh_cell(mesh, mesh_facet.entities(mesh.topology().dim())[0]);	  
	  // Get local index of facet with respect to the cell                 
	  uint local_facet = mesh_cell.index(mesh_facet);
	  
	  // Compute the lenght/area of the cell
	  real area = mesh_cell.volume();

	  // Compute a normal to the boundary
	  normal[0] += area * mesh_cell.normal(local_facet, 0);
	  normal[1] += area * mesh_cell.normal(local_facet, 1);
	  if (nsdim == 3)
	    normal[2] += area * mesh_cell.normal(local_facet, 2);

	  if(mesh.distdata().is_shared(id, 0)){
	    send_area.push_back(area);	    
	    for (int j = 0; j < Ncol; j++) 
	      send_buffn.push_back( mesh_cell.normal(local_facet, j) );
	  }
	}

      if(mesh.distdata().is_shared(id, 0)) {
	send_buff_indices.push_back(mesh.distdata().get_global(id, 0));
	send_buff.push_back(normal[0]);
	send_buff.push_back(normal[1]);
	if(nsdim == 3)
	  send_buff.push_back(normal[2]);
	
	send_ni.push_back( Nrow );
	send_ni.push_back( _noffset );
	send_ni.push_back( area_offset );
	
	// Init datastructures for shared data
	shared_noffset[mesh.distdata().get_global(id, 0)] = _offset;
	for(uint k=0; k< nsdim; k++)
	  shared_normal.push_back(0.0);
	_offset += nsdim;

	_noffset += Ncol * Nrow;
	area_offset += Nrow;
	num_cells[mesh.distdata().get_global(id, 0)] = Nrow;
      }      
    }

  }

  // Exchange values
  MPI_Status status;
  uint src,dest;
  int sh_count = send_buff_indices.size();
  int sh_count_data = send_buff.size();
  int sh_count_ndata = send_buffn.size();
  int sh_count_area = send_area.size();
  int recv_size_sh, recv_size_data, 
    recv_size_ndata, recv_count,recv_size_area; // recv_size_ni

  MPI_Barrier(dolfin::MPI::DOLFIN_COMM);
  MPI_Allreduce(&sh_count, &recv_size_sh, 1,MPI_INT,MPI_MAX, dolfin::MPI::DOLFIN_COMM);
  MPI_Allreduce(&sh_count_data, &recv_size_data, 
		1,MPI_INT,MPI_MAX, dolfin::MPI::DOLFIN_COMM);

  MPI_Allreduce(&sh_count_area, &recv_size_area, 
		1,MPI_INT,MPI_MAX, dolfin::MPI::DOLFIN_COMM);

  MPI_Allreduce(&sh_count_ndata, &recv_size_ndata, 
		1,MPI_INT,MPI_MAX, dolfin::MPI::DOLFIN_COMM);
  uint *recv_shared = new uint[ recv_size_sh];
  real *recv_data = new real[ recv_size_data ];
  real *recv_ndata = new real[ recv_size_ndata ];
  real *recv_area = new real[ recv_size_area];
  uint *recv_nidx = new uint[ 3* recv_size_sh ];

  for(int j=1; j < pe_size; j++){
    
    src = (rank - j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;
    
    MPI_Sendrecv(&send_buff_indices[0], sh_count, MPI_UNSIGNED, src, 1, 
		 recv_shared, recv_size_sh, MPI_UNSIGNED, dest, 1,
		 dolfin::MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status,MPI_UNSIGNED,&recv_count);

    MPI_Sendrecv(&send_buff[0], sh_count_data, MPI_DOUBLE, src, 1, 
		 recv_data, recv_size_data, MPI_DOUBLE, dest, 1,
		 dolfin::MPI::DOLFIN_COMM, &status);

    MPI_Sendrecv(&send_buffn[0], sh_count_ndata, MPI_DOUBLE, src, 1, 
		 recv_ndata, recv_size_ndata, MPI_DOUBLE, dest, 1,
		 dolfin::MPI::DOLFIN_COMM, &status);

    MPI_Sendrecv(&send_area[0], sh_count_area, MPI_DOUBLE, src, 1, 
		 recv_area, recv_size_area, MPI_DOUBLE, dest, 1,
		 dolfin::MPI::DOLFIN_COMM, &status);

    MPI_Sendrecv(&send_ni[0], (3*sh_count), MPI_UNSIGNED, src, 1, 
		 recv_nidx, 3*recv_size_sh, MPI_UNSIGNED, dest, 1,
		 dolfin::MPI::DOLFIN_COMM, &status);
    
    uint data_i = 0;
    uint data_ni = 0;
    for(int i=0; i<recv_count; i++) {
      uint glb_index = recv_shared[i];
      if(mesh.distdata().have_global(glb_index, 0) &&
	 glb_boundary.count(glb_index) > 0 &&
	 glb_boundary[glb_index]) {

	shared_normal[shared_noffset[glb_index]] += recv_data[data_i ]; 
	shared_normal[shared_noffset[glb_index]+1] += recv_data[data_i+1]; 
	if(nsdim == 3)
	  shared_normal[shared_noffset[glb_index]+2] += recv_data[data_i+2]; 
	
	num_cells[glb_index] += recv_nidx[ data_ni ];
	uint noffset = recv_nidx[ data_ni + 1];	
	for(uint k = 0; k < 3 * recv_nidx[data_ni]; k++) 
	  normal_block[ glb_index ].push_back( recv_ndata[ noffset++ ] );

	uint aoffset = recv_nidx[ data_ni + 2];
	for(uint k = 0; k < recv_nidx[ data_ni ]; k++)
	  shared_area_block[glb_index].push_back( recv_area[ aoffset++ ]);
      }
      
      // Increase data offsets
      data_i += nsdim;
      data_ni += 3;
      
    }
    

  }
  delete[] recv_shared;
  delete[] recv_data;
  delete[] recv_ndata;
  delete[] recv_nidx;  
  delete[] recv_area;
}
