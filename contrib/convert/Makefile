CFLAGS  = `pkg-config --cflags libxml-2.0` -I./  
LIBS    = `pkg-config --libs libxml-2.0` 
CXX     = cc

TARGET = convert
SRC = convert.c

all: $(TARGET)

clean:
	-rm -f *.o core *.core main.o $(OBJS) $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CFLAGS) $(SRC) $(LIBS) -o $@

