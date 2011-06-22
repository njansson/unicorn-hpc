CFLAGS  = `pkg-config --cflags unicorn` 
LIBS    = `pkg-config --libs unicorn` 
CXX     = `pkg-config --variable=compiler dolfin`

FORM_SRC = NSEContinuity3D.cpp NSEDualContinuity3D.cpp NSEDualGradient3D.cpp \
	NSEDualMomentum3D.cpp NSEMomentum3D.cpp NSEResContinuity3D.cpp NSEResMomentum3D.cpp Drag3D.cpp NSEMomentum2D.cpp NSEContinuity2D.cpp NSEDensity2D.cpp NSEDensity3D.cpp NavierStokesStress2D.cpp NavierStokesStress3D.cpp

FORM_OBJ = ${FORM_SRC:.cpp=.o}


all: $(FORM_OBJ)

install:

clean:
	-rm -f *.o core *.core $(FORM_OBJ) 
.cpp.o:
	$(CXX) $(CFLAGS) -c $< 
