SUBDIRS = unicorn
AM_CPPFLAGS = -I@srcdir@/unicorn
nobase_include_HEADERS = unicorn/ErrorEstimate.h  unicorn/init.h  unicorn/NodeNormal.h\
	  unicorn/SlipBC.h unicorn/NewSlipBC.h unicorn/NewDirichletBC.h \
	unicorn/SpaceTimeFunction.h  unicorn/TimeDependentPDE.h unicorn/UniParameters.h \
	unicorn/Project.h unicorn/EquiAffineMap.h unicorn/MeshQuality.h \
	unicorn/util.h unicorn/DataSetSampler.h unicorn/ElasticSmoother.h
lib_LIBRARIES = libunicorn.a
libunicorn_a_SOURCES = UniParameters.cpp NodeNormal.cpp SpaceTimeFunction.cpp\
	ErrorEstimate.cpp SlipBC.cpp init.cpp util.cpp \
	  NewSlipBC.cpp NewDirichletBC.cpp Project.cpp \
	TimeDependentPDE.cpp EquiAffineMap.cpp \
	MeshQuality.cpp DataSetSampler.cpp \
	Elasticity2D.cpp  ElasticityJac2D.cpp  ElasticityStress2D.cpp \
	Elasticity3D.cpp  ElasticityJac3D.cpp  ElasticityStress3D.cpp \
	ElasticSmoother.cpp


