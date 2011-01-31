// Copyright (C) 2010 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2010-09-13
// Last changed: 2011-01-18

#include <algorithm>
#include <limits>
#include <fstream>
#include <errno.h>
#include <sys/stat.h>
#include <dolfin.h>
#include <dolfin/io/BinaryFile.h>
#include <dolfin/fem/UFC.h>
#include "unicorn/AdaptiveRefinement.h"
#include "unicorn/AdaptiveRefinementProjectScalar.h"
#include "unicorn/AdaptiveRefinementProjectVector.h"

using namespace dolfin;
using namespace dolfin::unicorn;

//-----------------------------------------------------------------------------
void AdaptiveRefinement::refine(Mesh& mesh, MeshFunction<bool>& cell_marker)
				
{
  if(MPI::processNumber() == 0)
    dolfin_set("output destination","terminal");
  message("Adaptive refinement");
  message("cells before: %d", 
	  (dolfin::MPI::numProcesses() > 1 ? mesh.distdata().global_numCells() : mesh.numCells()));
  message("vertices before: %d", 
	  (dolfin::MPI::numProcesses() > 1 ? mesh.distdata().global_numVertices() : mesh.numVertices()));
  dolfin_set("output destination","silent");  
  
  MeshFunction<real> cell_refinement_marker_r(mesh);
  cell_refinement_marker_r.init(mesh.topology().dim());
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    cell_refinement_marker_r.set(c->index(), cell_marker(*c));
  }

  File refinefile("marked.pvd");
  refinefile << cell_refinement_marker_r;

  if(MPI::processNumber() == 0)
    dolfin_set("output destination","terminal");

  const std::string refine_type = dolfin_get("adapt_algorithm");
  if(refine_type == "rivara")
    RivaraRefinement::refine(mesh, cell_marker);
  else if(refine_type == "simple")
    mesh.refine(cell_marker, true);
  else
    dolfin::error("Unknown refinement algorithm");
  dolfin_set("output destination","silent");

  if(MPI::processNumber() == 0)
    dolfin_set("output destination","terminal");
  message("cells after: %d", 
	  (dolfin::MPI::numProcesses() > 1 ? mesh.distdata().global_numCells() : mesh.numCells()));
  message("vertices after: %d", 
	  (dolfin::MPI::numProcesses() > 1 ? mesh.distdata().global_numVertices() : mesh.numVertices()));
  dolfin_set("output destination","silent"); 
}
//-----------------------------------------------------------------------------
void AdaptiveRefinement::refine_and_project(Mesh& mesh, 			
					    std::vector<project_func> pf,	
					    MeshFunction<bool>& cell_marker)
{

  dolfin_set("Load balancer redistribute", false);

  if(MPI::processNumber() == 0)
    dolfin_set("output destination","terminal");
  message("Adaptive refinement (with projection)");
  message("cells before: %d", mesh.distdata().global_numCells());
  message("vertices before: %d", mesh.distdata().global_numVertices());
  
  const std::string refine_type = dolfin_get("adapt_algorithm");
  if(refine_type == "simple")
    LoadBalancer::balance(mesh, cell_marker);
  else if(refine_type == "rivara")
    LoadBalancer::balance(mesh, cell_marker, LoadBalancer::LEPP);
  else
    dolfin::error("Unknown refinement algorithm");
  //  dolfin_set("output destination","silent");

  MeshFunction<real> cell_refinement_marker_r(mesh);
  cell_refinement_marker_r.init(mesh.topology().dim());
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    cell_refinement_marker_r.set(c->index(), cell_marker(*c));
  }
  
  File refinefile("marked.pvd");
  refinefile << cell_refinement_marker_r;

  MeshFunction<uint> *partitions = mesh.data().meshFunction("partitions");
   
  real *x_values, *y_values, *z_values;
  uint *x_rows, *y_rows, *z_rows;
  uint x_m, y_m, z_m;
  x_values = y_values = z_values = 0;
  x_rows = y_rows, z_rows = 0;
  x_m = y_m = z_m = 0;
  Function coarse_x, coarse_y, coarse_z;
  Vector x_coarse_x, x_coarse_y, x_coarse_z;
  Form *pre_x = new AdaptiveRefinementProjectScalarLinearForm(coarse_x);
  Form *pre_y = new AdaptiveRefinementProjectScalarLinearForm(coarse_y);
  Form *pre_z = new AdaptiveRefinementProjectScalarLinearForm(coarse_z);
  coarse_x.init(mesh, x_coarse_x, *pre_x, 1);
  coarse_y.init(mesh, x_coarse_y, *pre_y, 1);
  coarse_z.init(mesh, x_coarse_z, *pre_z, 1);

  
  for (std::vector<project_func>::iterator it = pf.begin(); 
       it != pf.end(); ++it)
  {
    if(it->first->type() != Function::discrete)
      error("Projection only implemented for discrete functions");   

    AdaptiveRefinement::decompose_func(mesh, it->first, 
				       it->second.second, *(it->second.first), 
				       coarse_x, coarse_y, coarse_z);

    AdaptiveRefinement::redistribute_func(mesh, &coarse_x, &x_values, 
					  &x_rows, x_m, *partitions);
    AdaptiveRefinement::redistribute_func(mesh, &coarse_y, &y_values, 
					  &y_rows, y_m, *partitions);
    
    AdaptiveRefinement::redistribute_func(mesh, &coarse_z, &z_values, 
					  &z_rows, z_m, *partitions);
        
  }

  MeshFunction<bool> new_cell_marker;
  mesh.distribute(*partitions, cell_marker, new_cell_marker);

  Mesh new_mesh;
  new_mesh = mesh;
  RivaraRefinement::refine(new_mesh, new_cell_marker, 0.0, 0.0, 0.0, false);
  new_mesh.renumber();

  
  if(MPI::processNumber() == 0)
    if(mkdir("../scratch", S_IRWXU) < 0)
      perror("mkdir failed");

  int p_count = 0;
  for (std::vector<project_func>::iterator it = pf.begin(); 
       it != pf.end(); ++it)
  {
    if(it->first->type() != Function::discrete)
      error("Projection only implemented for discrete functions");   

    Function post_x, post_y, post_z;
    Vector x_post_x, x_post_y, x_post_z;
    post_x.init(mesh, x_post_x, *pre_x, 1);
    post_y.init(mesh, x_post_y, *pre_y, 1);
    post_z.init(mesh, x_post_z, *pre_z, 1);

    post_x.vector().set(x_values, x_m, x_rows);
    post_y.vector().set(y_values, y_m, y_rows);
    post_z.vector().set(z_values, z_m, z_rows);

    post_x.sync_ghosts();
    post_y.sync_ghosts();
    post_z.sync_ghosts();

    delete[] x_values;
    delete[] y_values;
    delete[] z_values;

    delete[] x_rows;
    delete[] y_rows;
    delete[] z_rows;



    Vector xproj;    
    AdaptiveRefinement::project(new_mesh, post_x, post_y, post_z, xproj);

    std::stringstream p_filename;
    p_filename << "../scratch/projected_" << p_count++ << "_" << MPI::processNumber() << ".bin" << std::ends;
    File p_file(p_filename.str());
    p_file << xproj;


    delete pre_x, pre_y, pre_z;

  }

  mesh = new_mesh;  
  mesh.renumber();

}
//-----------------------------------------------------------------------------
void AdaptiveRefinement::redistribute_func(Mesh& mesh, Function *f, 
					   real **vp, uint **rp, uint& m,
					   MeshFunction<uint>& distribution)
{


  uint pe_rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();
  uint target_proc, src, dest, recv_size, local_size;
  
  MPI_Status status;

  real *values = new real[f->vector().local_size()];
  f->vector().get(values);
  
  Array<real> *send_buffer = new Array<real>[pe_size];
  Array<uint> *send_buffer_indices = new Array<uint>[pe_size];


  MeshFunction<bool> marked(mesh, 0);
  marked = false;

  std::vector<std::pair<uint, real> > recv_data;
  local_size = 0;

  uint nsdim = mesh.topology().dim();
  real value;
  uint global_index;

  for (CellIterator c(mesh); !c.end(); ++c) 
  {

    target_proc = distribution.get(*c);

    for (VertexIterator v(*c); !v.end(); ++v)
    {
      
      global_index = mesh.distdata().get_global(v->index(), 0);
      f->vector().get(&value, 1, &global_index);
      
      if (target_proc == pe_rank && 
	  !mesh.distdata().is_ghost(v->index(), 0) &&
	  !marked.get(*v))
      {	
	
	std::pair<uint, real> p(global_index, value);
	recv_data.push_back(p);	
	marked.set(*v, true);
	continue;	
      }
      
      if (!mesh.distdata().is_ghost(v->index(), 0) && !marked.get(*v))
      {	
	
	send_buffer[target_proc].push_back(value);
	send_buffer_indices[target_proc].push_back(global_index);	
	marked.set(*v, true);
      }
    }
    local_size = std::max(local_size, (uint) send_buffer[target_proc].size());
  }

  delete[] values;

  MPI_Allreduce(&local_size, &recv_size, 1, 
		MPI_UNSIGNED, MPI_MAX, MPI::DOLFIN_COMM);

  real *recv_buffer = new real[recv_size];
  uint *recv_buffer_indices = new uint[recv_size];


  
  int recv_count = 0;

  for (uint j = 1; j < pe_size; j++)
  {
    src = (pe_rank - j + pe_size) % pe_size;
    dest = (pe_rank + j) % pe_size;

    MPI_Sendrecv(&send_buffer[dest][0], send_buffer[dest].size(), MPI_DOUBLE,
		 dest, 1, recv_buffer, recv_size, MPI_DOUBLE, src, 1,
		 MPI::DOLFIN_COMM, &status);

    MPI_Sendrecv(&send_buffer_indices[dest][0], send_buffer_indices[dest].size(), 
		 MPI_UNSIGNED, dest, 1, recv_buffer_indices, recv_size, 
		 MPI_UNSIGNED, src, 1, MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

      
    for (int i = 0; i < recv_count; i++) 
    {
      std::pair<uint, real> p(recv_buffer_indices[i], recv_buffer[i]);
      recv_data.push_back(p);
    }    

  }

  delete[] recv_buffer;
  delete[] recv_buffer_indices;

  for (uint i = 0; i < pe_size; i++) 
  {
    send_buffer[i].clear();
    send_buffer_indices[i].clear();
  }
  
  delete[] send_buffer;
  delete[] send_buffer_indices;

  local_size = recv_data.size();
  values = new real[local_size];
  uint *rows = new uint[local_size];

  for (uint i = 0; i < local_size; i++) {
    rows[i] = recv_data[i].first;
    values[i] = recv_data[i].second;    
  }

  *vp = values;
  *rp = rows;
  m = local_size;
  
}
//-----------------------------------------------------------------------------
void AdaptiveRefinement::decompose_func(Mesh& mesh, Function *f, uint offset, 
					Form& form, Function &f_x, 
					Function &f_y, Function &f_z)
{

  UFC ufc(form.form(), mesh, form.dofMaps());  
  uint local_dim = (form.dofMaps())[offset].local_dimension();
  uint *indices = new uint[local_dim];
  uint new_index;
  MeshFunction<bool> marked(mesh, 0);
  marked = false;
  
  real dof_value;
  for (CellIterator c(mesh); !c.end(); ++c) 
  {
    
    ufc.update(*c, mesh.distdata());
    (form.dofMaps())[offset].tabulate_dofs(indices, ufc.cell, c->index());
    
    for (VertexIterator v(*c); !v.end(); ++v)
    {
      
      uint *cvi = c->entities(0);
      uint ci = 0;
      for(ci = 0; ci < c->numEntities(0); ci++)
	if(cvi[ci] == v->index())
	  break;
      
      if (!mesh.distdata().is_ghost(v->index(), 0) && !marked.get(*v))
      {

	new_index = mesh.distdata().get_global(v->index(), 0);

	f->vector().get(&dof_value, 1, &indices[ci]);
	f_x.vector().set(&dof_value, 1, &new_index);

	f->vector().get(&dof_value, 1, &indices[ci + c->numEntities(0)]);
	f_y.vector().set(&dof_value, 1, &new_index);

	f->vector().get(&dof_value, 1, &indices[ci + 2 * c->numEntities(0)]);
	f_z.vector().set(&dof_value, 1, &new_index);
	
	marked.set(*v, true);
	continue;
	
      }
    }
  }

  f_x.sync_ghosts();
  f_y.sync_ghosts();
  f_z.sync_ghosts();


  delete[] indices;

}
//-----------------------------------------------------------------------------
void AdaptiveRefinement::project(Mesh& new_mesh, Function& post_x,
				 Function& post_y, Function& post_z, 
				 Vector& x_proj)
{
  
  Function projected;
  Form *refined = new AdaptiveRefinementProjectVectorLinearForm(projected);
  projected.init(new_mesh, x_proj, *refined, 0);
  
  
  UFC ufc(refined->form(), new_mesh, refined->dofMaps());
  
  real test_value;
  real x[3];
  real *vv = new real[x_proj.local_size()];
  uint *indices = new uint[x_proj.local_size()];
  uint *local_indices = new uint[refined->dofMaps()[0].local_dimension()];
  uint i = 0;
  MeshFunction<bool> processed(new_mesh, 0);
  processed = false;
  
  projected.vector().zero();
  projected.sync_ghosts();
  
  for (CellIterator c(new_mesh); !c.end(); ++c) {
    
    
    ufc.update(*c, new_mesh.distdata());
    (refined->dofMaps())[0].tabulate_dofs(local_indices, ufc.cell, c->index());
    
    for (VertexIterator v(*c); !v.end(); ++v) {
      
      
      uint *cvi = c->entities(0);
      uint ci = 0;
      for(ci = 0; ci < c->numEntities(0); ci++)
	if(cvi[ci] == v->index())
	  break;
      
      if(new_mesh.distdata().is_ghost(v->index(), 0) || processed.get(*v))
	continue;
      processed.set(*v, true);
      
      Vertex *v_e = 0;

      x[0] = v->x()[0];
      x[1] = v->x()[1];
      x[2] = v->x()[2];           
      post_x.eval(&test_value, &x[0]);      
      if (test_value == std::numeric_limits<real>::infinity()) {
	for (EdgeIterator e(*v); !e.end(); ++e) {
	  const uint *edge_v = e->entities(0);

	  (edge_v[0] != v->index() ? 
	   v_e = new Vertex(new_mesh, edge_v[0]) : v_e = new Vertex(new_mesh, edge_v[1]));
	  x[0] = v_e->x()[0];
	  x[1] = v_e->x()[1];
	  x[2] = v_e->x()[2];        
	  post_x.eval(&test_value, &x[0]);   
	if (test_value != std::numeric_limits<real>::infinity()) 
	    break;
	}

	if (test_value == std::numeric_limits<real>::infinity()) 
	  error("Couldn't find any suitable projection point");	   	
      }

      vv[i] = test_value;
      indices[i++] = local_indices[ci];
      
      
      post_y.eval(&test_value, &x[0]);      
      if (test_value != std::numeric_limits<real>::infinity())  {
	vv[i] = test_value;
	indices[i++] = local_indices[ci +  c->numEntities(0)];
      }
      
      post_z.eval(&test_value, &x[0]);      
      if (test_value != std::numeric_limits<real>::infinity())  {
	vv[i] = test_value;
	indices[i++] = local_indices[ci + 2 * c->numEntities(0)]; 
      }

      if (v_e)
	delete v_e;
    }
  }
  
  if ( i > 0) {
    x_proj.set(vv, i, indices);
    x_proj.apply();
  }
  projected.sync_ghosts();
  
  delete[] vv;
  delete[] indices;
  delete[] local_indices;
}
//-----------------------------------------------------------------------------
