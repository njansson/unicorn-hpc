CFLAGS  = `pkg-config --cflags unicorn` -I./ -I../
LIBS    = `pkg-config --libs unicorn` 
CXX     = `pkg-config --variable=compiler dolfin`

TARGET  = cylinder
FORM_SRC = ../unicorn/NSEContinuity3D.cpp ../unicorn/NSEDualContinuity3D.cpp \
	 ../unicorn/NSEDualGradient3D.cpp ../unicorn/NSEDualMomentum3D.cpp \
	 ../unicorn/NSEMomentum3D.cpp ../unicorn/NSEResContinuity3D.cpp  \
	../unicorn/NSEResMomentum3D.cpp ../unicorn/Drag3D.cpp \
	 ../unicorn/NSEMomentum2D.cpp ../unicorn/NSEContinuity2D.cpp
SRC	 = ../NSESolver.cpp

FORM_OBJ = ${FORM_SRC:.cpp=.o}
OBJECTS = ${SRC:.cpp=.o}

all: $(TARGET)

install:

clean:
	-rm -f *.o core *.core main.o main_slip.o $(FORM_OBJ) $(OBJECTS) $(TARGET)

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

