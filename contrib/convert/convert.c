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
#include <libxml/xmlreader.h>

int parse_header(xmlTextReaderPtr xml_reader) {

  printf("celltype=%s dim=%s\n",
	 xmlTextReaderGetAttribute(xml_reader,"celltype"),
	 xmlTextReaderGetAttribute(xml_reader,"dim"));

  return 0;
}


int parse_vertices(xmlTextReaderPtr xml_reader) {

  int i, size;
  double x,y,z;

  size = atoi((const char *)xmlTextReaderGetAttribute(xml_reader,"size"));

  x = y = z = 0.0;
  xmlTextReaderRead(xml_reader);    
  xmlTextReaderRead(xml_reader);    
  for (i = 0 ; i < size;
       xmlTextReaderRead(xml_reader),
	 xmlTextReaderRead(xml_reader), i++) {
    
    x = atof((const char *)xmlTextReaderGetAttribute(xml_reader,"x"));
    y = atof((const char *)xmlTextReaderGetAttribute(xml_reader,"y"));
    if ( size == 3)
      z = atof((const char *)xmlTextReaderGetAttribute(xml_reader,"z"));
    
    printf("x:%g y:%g z:%g\n", x,y,z);
        
  }  
  
  return 0;
}

int parse_cells(xmlTextReaderPtr xml_reader) {

  int i, size;
  int v0, v1, v2, v3;
  

  size = atoi((const char *)xmlTextReaderGetAttribute(xml_reader,"size"));
  
  xmlTextReaderRead(xml_reader);    
  xmlTextReaderRead(xml_reader);    
  for (i = 0 ; i < size;
       xmlTextReaderRead(xml_reader),
	 xmlTextReaderRead(xml_reader), i++) {
    
    v0 = atoi((const char *)xmlTextReaderGetAttribute(xml_reader,"v0"));
    v1 = atoi((const char *)xmlTextReaderGetAttribute(xml_reader,"v1"));
    v2 = atoi((const char *)xmlTextReaderGetAttribute(xml_reader,"v2"));
    
    printf("v0:%d v1:%d v2:%d\n", v0,v1,v2);
        
  }  
  
  return 0;
}




int main(int argc, char *argv[]) {
  
  FILE  *binary_fp;
  xmlTextReaderPtr xml_reader;
  int status;

  if (argc < 2 ) { 
    fprintf(stderr, "Usage: ./convert <xml mesh> <binary mesh>\n");
    return -1;
  }
  
  xml_reader = xmlNewTextReaderFilename(argv[1]);
  
  if(xml_reader == NULL) {
    fprintf(stderr, "Cant open DOLFIN xml file %s\n", argv[1]);
    return -1;
  }
  
  /* Parse header */
  status = xmlTextReaderRead(xml_reader);
  status = xmlTextReaderRead(xml_reader);
  status = xmlTextReaderRead(xml_reader);
  parse_header(xml_reader);
  
  /* Parse vertices */
  status = xmlTextReaderRead(xml_reader);
  status = xmlTextReaderRead(xml_reader);
  parse_vertices(xml_reader);

  /* Parse vertices */
  status = xmlTextReaderRead(xml_reader);
  status = xmlTextReaderRead(xml_reader);
  parse_cells(xml_reader);


    
  xmlFreeTextReader(xml_reader);


  //  binary_fp = fopen(argv[2], "w");
 

   

  


  //  fclose(xml_fp);
  //  fclose(binary_fp);

 
  return 0;
}
