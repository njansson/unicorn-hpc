/*
 * ----------------------------------------------------------------------------
 * "THE BEER-WARE LICENSE" (Revision 42):
 * <njansson@csc.kth.se> wrote this file. As long as you retain this notice you
 * can do whatever you want with this stuff. If we meet some day, and you think
 * this stuff is worth it, you can buy me a beer in return Niclas Jansson
 * ----------------------------------------------------------------------------
 */

#include <cstdlib>
#include <fstream>
#include <list>
#include <map>
#include <string>
#include <vector>
#include <dolfin.h>

using namespace dolfin;


struct dVertex {
  dVertex(int _id, Point _p) 
  { id = _id; p = _p;};
  unsigned int id;
  Point p;
};

struct dCell {
  dCell(Array<unsigned int> _vertices) 
  { vertices = _vertices;};
  Array<unsigned int> vertices;
};

int main(int argc, char *argv[]) {
  
  if (argc < 3 ) {
    std::cout<< "Usage: ./extract <checkpoint> <num>" << std::endl;
    return 0;
  }
  
  std::ostringstream _fname;
  _fname << argv[1] << "_";

  std::vector<std::string> chkp_files;

  for (int i = 0; i <  atoi(argv[2]); i++) {
    std::ostringstream tmp;
    tmp << _fname.str();
    tmp << i << ".chkp";
    chkp_files.push_back(tmp.str());
  }
  

  CellType::Type type;
  double t;
  unsigned int id, tdim, gdim, num_vertices, num_cells, num_entities;
  std::ifstream in;
  std::list<dVertex> vlist;
  std::list<dCell> dmesh;
  for(std::vector<std::string>::iterator it = chkp_files.begin(); 
      it != chkp_files.end(); ++it) { 

    in.open(it->c_str(), std::ifstream::binary);

    if(in.good()) {
      in.read((char *)&id, sizeof(unsigned int));
      in.read((char *)&t, sizeof(double));
      
      in.read((char *)&type, sizeof(CellType::Type));
      in.read((char *)&tdim, sizeof(unsigned int));
      in.read((char *)&gdim, sizeof(unsigned int));

      in.read((char *)&num_vertices, sizeof(unsigned int));
      in.read((char *)&num_cells, sizeof(unsigned int));
      in.read((char *)&num_entities, sizeof(unsigned int));

      double *coords = new double[gdim *num_vertices];
      in.read((char *)coords, (gdim * num_vertices) * sizeof(double));
      
      std::list<dVertex> _vlist;
      int vi = 0;
      for(unsigned int i = 0; i < gdim * num_vertices; i += gdim) {
	switch(gdim)
	{
	case 2:
	  _vlist.push_back(dVertex(vi++,  Point(coords[i], coords[i+1]))); break;
	case 3:
	  _vlist.push_back(dVertex(vi++,  Point(coords[i], coords[i+1], coords[i+2]))); break;
	}
      }
      delete[] coords;

      unsigned int *cells = new unsigned int[num_entities * num_cells];
      in.read((char *)cells, (num_entities * num_cells) * sizeof(unsigned int));

      unsigned int *mapping = new unsigned int[_vlist.size()];
      in.read((char *)mapping, _vlist.size() * sizeof(unsigned int));

      unsigned int *mp = &mapping[0];
      std::map<unsigned int, unsigned int> vmap;      
      for(std::list<dVertex>::iterator it = _vlist.begin(); 
	  it != _vlist.end(); ++it) 
	vmap[it->id] = *(mp++);
      delete[] mapping;	
      
      Array<unsigned int> v;
      for(unsigned int i = 0; i < num_entities * num_cells; i += num_entities) {
	v.clear();
	for(unsigned int j = 0; j < num_entities; j++)
	  v.push_back(vmap[cells[i+j]]);
	dmesh.push_back(dCell(v));
      }
      delete[] cells;

      unsigned int num_ghost;
      in.read((char *)&num_ghost, sizeof(unsigned int));
      unsigned int *ghosts = new unsigned int[2 * num_ghost];
      in.read((char *)ghosts, 2*num_ghost * sizeof(unsigned int));
      
      std::set<unsigned int> ghost_set;
      for (unsigned int i = 0; i < 2 * num_ghost; i += 2)
	ghost_set.insert(ghosts[i]);
      delete[] ghosts;

      for(std::list<dVertex>::iterator it = _vlist.begin(); 
	  it != _vlist.end(); ++it) 
	if(ghost_set.find(it->id) == ghost_set.end()) 
	  vlist.push_back(dVertex(vmap[it->id], it->p));
    }
    in.close();
  }

  Mesh mesh;
  MeshEditor editor;
  editor.open(mesh, type, tdim, gdim);
  editor.initVertices(vlist.size());
  editor.initCells(dmesh.size());

  
  for(std::list<dVertex>::iterator it = vlist.begin(); it != vlist.end(); ++it)
    editor.addVertex(it->id, it->p);

  unsigned int ci = 0;
  for(std::list<dCell>::iterator it = dmesh.begin(); it != dmesh.end(); ++it)
    editor.addCell(ci++, it->vertices);
  
  editor.close();

  File mesh_file("mesh.xml");
  mesh_file << mesh;

}
