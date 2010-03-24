CFLAGS  = `pkg-config --cflags unicorn` 
CXX     = `pkg-config --variable=compiler dolfin`

SRC = NSESolver.cpp

OBJ = ${SRC:.cpp=.o}


all: $(OBJ)

clean:
	-rm -f *.o core *.core $(FORM_OBJ) 
.cpp.o:
	$(CXX) $(CFLAGS) -c $< 
