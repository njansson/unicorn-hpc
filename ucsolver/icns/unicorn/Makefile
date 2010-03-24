CFLAGS  = `pkg-config --cflags unicorn` 
CXX     = `pkg-config --variable=compiler dolfin`

FORM_SRC = NSEContinuity3D.cpp NSEDualContinuity3D.cpp NSEDualGradient3D.cpp \
	NSEDualMomentum3D.cpp NSEMomentum3D.cpp NSEResContinuity3D.cpp NSEResMomentum3D.cpp \
	Drag3D.cpp

FORM_OBJ = ${FORM_SRC:.cpp=.o}


all: $(FORM_OBJ)

install:

clean:
	-rm -f *.o core *.core $(FORM_OBJ) 
.cpp.o:
	$(CXX) $(CFLAGS) -c $< 
