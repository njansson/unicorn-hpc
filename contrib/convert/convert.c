/*
 * ----------------------------------------------------------------------------
 * "THE BEER-WARE LICENSE" (Revision 42):
 * <njansson@csc.kth.se> wrote this file. As long as you retain this notice you
 * can do whatever you want with this stuff. If we meet some day, and you think
 * this stuff is worth it, you can buy me a beer in return Niclas Jansson
 * ----------------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libxml/xmlreader.h>

void progress(int *state) {
  switch(*state)
    {
    case 0:
      putchar(45);
      break;
    case 1:      
      putchar(92);
      break;
    case 2:
      putchar(124);
      break;
    case 3:
      putchar(47);
      break;
    }
  fflush(stdout);
  putchar(8);
  *state = (*state + 1)%4;
}
int parse_header(xmlTextReaderPtr xml_reader, FILE *binary_fp,
		 int *dim, int *celltype) {

  *dim = atoi((const char *)xmlTextReaderGetAttribute(xml_reader,"dim"));
	 
  if(strcmp((const char *) xmlTextReaderGetAttribute(xml_reader,"celltype"),
	     "triangle") == 0)
    *celltype = 0;
  else if (strcmp((const char *)
		   xmlTextReaderGetAttribute(xml_reader,"celltype"),
		   "tetrahedron") == 0)
    *celltype = 1;

  fwrite(dim, sizeof(int), 1, binary_fp);
  fwrite(celltype, sizeof(int), 1, binary_fp);
  return 0;
}


int parse_vertices(xmlTextReaderPtr xml_reader, FILE *binary_fp, int dim) {

  int i, size,state;
  double *data, *dp;

  size = atoi((const char *)xmlTextReaderGetAttribute(xml_reader,"size"));
  data = malloc(size * dim * sizeof(double));
  dp = &data[0];

  state = 0;
  printf("Reading vertices...");
  xmlTextReaderRead(xml_reader);    
  xmlTextReaderRead(xml_reader);    
  for (i = 0 ; i < size;
       xmlTextReaderRead(xml_reader),
	 xmlTextReaderRead(xml_reader), i++) {
    progress(&state);
    *(dp++) = atof((const char *)xmlTextReaderGetAttribute(xml_reader,"x"));
    *(dp++) = atof((const char *)xmlTextReaderGetAttribute(xml_reader,"y"));
    if ( dim == 3)
      *(dp++) = atof((const char *)xmlTextReaderGetAttribute(xml_reader,"z"));
    
  }  
  printf("Done\n");
  fwrite(&size, sizeof(int), 1, binary_fp);
  fwrite(data, sizeof(double), size * dim, binary_fp);
  free(data);
  return 0;
}

int parse_cells(xmlTextReaderPtr xml_reader, FILE *binary_fp, int celltype) {

  int i, size, state;
  int *data, *dp;

  size = atoi((const char *)xmlTextReaderGetAttribute(xml_reader,"size"));
  data = malloc(size * (3 + celltype) * sizeof(int));
  dp = &data[0];

  state = 0;
  printf("Reading cells...");  
  xmlTextReaderRead(xml_reader);    
  xmlTextReaderRead(xml_reader);    
  for (i = 0 ; i < size;
       xmlTextReaderRead(xml_reader),
	 xmlTextReaderRead(xml_reader), i++) {
    progress(&state);
    *(dp++) = atoi((const char *)xmlTextReaderGetAttribute(xml_reader,"v0"));
    *(dp++) = atoi((const char *)xmlTextReaderGetAttribute(xml_reader,"v1"));
    *(dp++) = atoi((const char *)xmlTextReaderGetAttribute(xml_reader,"v2"));
    if (celltype == 1)
      *(dp++) = atoi((const char *)xmlTextReaderGetAttribute(xml_reader,"v3"));           
  }  
  printf("Done\n");
  fwrite(&size, sizeof(int), 1, binary_fp);
  fwrite(data, sizeof(int), size * (3 + celltype), binary_fp);
  free(data);

  return 0;
}


int main(int argc, char *argv[]) {
  
  FILE  *binary_fp;
  xmlTextReaderPtr xml_reader;
  int dim, celltype;

  if (argc < 3 ) { 
    fprintf(stderr, "Usage: ./convert <xml mesh> <binary mesh>\n");
    return -1;
  }
  
  xml_reader = xmlNewTextReaderFilename(argv[1]);
  
  if(xml_reader == NULL) {
    fprintf(stderr, "Cant open DOLFIN xml file %s\n", argv[1]);
    return -1;
  }

  binary_fp = fopen(argv[2], "w");
  
  printf("Converting DOLFIN-xml to flat binary\n");

  /* Parse header */
  xmlTextReaderRead(xml_reader);
  xmlTextReaderRead(xml_reader);
  xmlTextReaderRead(xml_reader);
  parse_header(xml_reader, binary_fp, &dim, &celltype);
    
  /* Parse vertices */
  xmlTextReaderRead(xml_reader);
  xmlTextReaderRead(xml_reader);
  parse_vertices(xml_reader, binary_fp, dim);

  /* Parse cells */
  xmlTextReaderRead(xml_reader);
  xmlTextReaderRead(xml_reader);
  parse_cells(xml_reader, binary_fp, celltype);

  xmlFreeTextReader(xml_reader);
 
  fclose(binary_fp);
  return 0;
}
