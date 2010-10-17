// Copyright (C) 2010 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2010-09-13
// Last changed: 2010-10-11

#include <algorithm>
#include <limits>
#include <fstream>
#include <dolfin/io/BinaryFile.h>
#include <dolfin.h>
#include "unicorn/AdaptiveRefinement.h"

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
  message("Adaptive refinement");
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

  MeshFunction<uint> *partitions = mesh.data().meshFunction("partitions");
   
  real *values = 0;
  uint *rows = 0;
  uint m = 0;
  for (std::vector<project_func>::iterator it = pf.begin(); 
       it != pf.end(); ++it)
  {
    if(it->first->type() != Function::discrete)
      error("Projection only implemented for discrete functions");   
    //    (*it)->redistribute(*partitions);    
    AdaptiveRefinement::redistribute_func(mesh, it->first, 
					  &values, &rows, m, *partitions);
    if(values == 0)
      error("pre Kaos");
    File test("pre_func.pvd");
    test << *(it->first);
  }

  MeshFunction<bool> new_cell_marker;
  mesh.distribute(*partitions, cell_marker, new_cell_marker);
  //  mesh.renumber();

  Mesh new_mesh;
  new_mesh = mesh;
  RivaraRefinement::refine(new_mesh, cell_marker, 0.0,0.0,0.0, false);
  new_mesh.renumber();

  for (std::vector<project_func>::iterator it = pf.begin(); 
       it != pf.end(); ++it)
  {
    if(it->first->type() != Function::discrete)
      error("Projection only implemented for discrete functions");   
    Function tmp;
    Vector xtmp;
    message("offset %d \n\n\n\n\n\n\n", it->second.second);
    tmp.init(mesh, xtmp, *(it->second.first), it->second.second);
    message("local size %d\n\n\n", xtmp.local_size());
    //    tmp.vector().set(values);
    tmp.vector().set(values, m , rows);
    tmp.vector().apply();
    tmp.sync_ghosts();

    if(values == 0)
      error("Kaos");
    File test("post_func.pvd");
    test << tmp;
    //test << *(it->first);

    Vector xtmp_new;
    xtmp_new.init(new_mesh.numVertices() - new_mesh.distdata().num_ghost(0));

    real test_value;    
    real x[3];
    real *vv = new real[new_mesh.numVertices()];
    uint *indices = new uint[new_mesh.numVertices()];
    uint i = 0;
    for (VertexIterator v(new_mesh); !v.end(); ++v) {
      if(new_mesh.distdata().is_ghost(v->index(), 0))
	continue;
      x[0] = v->x()[0];
      x[1] = v->x()[1];
      x[2] = v->x()[2];
      tmp.eval(&test_value, &x[0]);
      if (!std::isinf(test_value)) {
	message("%g %g %g: %g",x[0], x[1], x[2] , test_value);
	vv[i] = test_value;
	indices[i++] = new_mesh.distdata().get_global(v->index(), 0);
      }

    }
    xtmp_new.set(vv, i, indices);
    xtmp_new.apply();

    BinaryFile bin_test("test.bin");
    bin_test << xtmp_new;
  }

  

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

  for (CellIterator c(mesh); !c.end(); ++c) 
  {

    target_proc = distribution.get(*c);
    

    for (VertexIterator v(mesh); !v.end(); ++v)
    {


      if (target_proc == pe_rank && 
	  !mesh.distdata().is_ghost(v->index(), 0) &&
	  !marked.get(*v))
      {
	std::pair<uint, real> p(mesh.distdata().get_global(v->index(),0),values[v->index()]);
	recv_data.push_back(p);
	marked.set(*v, true);
	continue;
      }

      if (!mesh.distdata().is_ghost(v->index(), 0) && !marked.get(*v))
      {
	send_buffer[target_proc].push_back(values[v->index()]);
	send_buffer_indices[target_proc].push_back(mesh.distdata().get_global(v->index(), 0));
	marked.set(*v, true);
      }
    }
    local_size = std::max(local_size, send_buffer[target_proc].size());
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


  less_pair comp;
  std::sort(recv_data.begin(), recv_data.end(), comp);

  local_size = recv_data.size();
  if(local_size == 0)
    error("Kaos redistribute");
  else
    message("realloc with size %d", local_size);
  values = new real[local_size];
  uint *rows = new uint[local_size];

  message("pre local_size %d", local_size);
  for (uint i = 0; i < local_size; i++)  {
    rows[i] = recv_data[i].first;
    values[i] = recv_data[i].second;    
  }

  *vp = values;
  *rp = rows;
  m = local_size;

  
  /*
  std::ofstream out("p.data", std::ofstream::binary);
  out.write((char *)&local_size, sizeof(uint));
  out.write((char *)values, local_size * sizeof(real));
  
  out.close();
  */


}
//-----------------------------------------------------------------------------
