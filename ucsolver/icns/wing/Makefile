CFLAGS  = `pkg-config --cflags unicorn` -I./ -I../
LIBS    = `pkg-config --libs unicorn` 
CXX     = `pkg-config --variable=compiler dolfin`

TARGET  = wing
FORM_SRC = ../unicorn/NSEContinuity3D.cpp ../unicorn/NSEDualContinuity3D.cpp \
	 ../unicorn/NSEDualGradient3D.cpp ../unicorn/NSEDualMomentum3D.cpp \
	 ../unicorn/NSEMomentum3D.cpp ../unicorn/NSEResContinuity3D.cpp  \
	../unicorn/NSEResMomentum3D.cpp ../unicorn/Drag3D.cpp	
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

$(TARGET): forms  solver
	$(CXX) $(CFLAGS) main.cpp $(FORM_OBJ) $(OBJECTS) $(LIBS) -o $@

.cpp.o:
	$(CXX) $(CFLAGS) -c $< 

