CFLAGS  = `pkg-config --cflags unicorn` -I./ -I../
LIBS    = `pkg-config --libs unicorn` 
CXX     = `pkg-config --variable=compiler dolfin`

TARGET  = elastic_test

all: $(TARGET)

install:

clean:
	-rm -f *.o core *.core main.o $(TARGET)

forms:	$(FORM_SRC)
	$(MAKE) -C ../unicorn

solver: $(SOLVER_SRC)
	$(MAKE) -C ../


wedge.xml:
	wget http://www.csc.kth.se/~jjan/private/wedge.xml

$(TARGET): main.cpp wedge.xml
	$(CXX) $(CFLAGS) -g main.cpp $(LIBS) -o $@

.cpp.o:
	$(CXX) $(CFLAGS) -c $< 

