CFLAGS  = `pkg-config --cflags dolfin` -I./  
LIBS    = `pkg-config --libs dolfin` 
CXX     = `pkg-config --variable=compiler dolfin`

TARGET = extract
SRC = extract.cpp

all: $(TARGET)

clean:
	-rm -f *.o core *.core main.o $(OBJS) $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CFLAGS) $(SRC) $(LIBS) -o $@

