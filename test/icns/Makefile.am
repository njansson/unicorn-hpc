AM_LDFLAGS = -L$(top_builddir)/src -lunicorn @DOLFIN_LIBS@ 
AM_CPPFLAGS = @DOLFIN_CFLAGS@ -I@top_builddir@/ucsolver/icns -I@top_builddir@/src
check_PROGRAMS = cylinder
cylinder_SOURCES = main.cpp \
	$(top_builddir)/ucsolver/icns/unicorn/NSEContinuity3D.cpp \
	$(top_builddir)/ucsolver/icns/unicorn/NSEDualContinuity3D.cpp\
	$(top_builddir)/ucsolver/icns/unicorn/NSEDualGradient3D.cpp \
        $(top_builddir)/ucsolver/icns/unicorn/NSEDualMomentum3D.cpp \
	$(top_builddir)/ucsolver/icns/unicorn/NSEMomentum3D.cpp \
	$(top_builddir)/ucsolver/icns/unicorn/NSEResContinuity3D.cpp \
	$(top_builddir)/ucsolver/icns/unicorn/NSEResMomentum3D.cpp \
        $(top_builddir)/ucsolver/icns/unicorn/Drag3D.cpp \
        $(top_builddir)/ucsolver/icns/unicorn/Drag3D_Int.cpp \
        $(top_builddir)/ucsolver/icns/NSESolver.cpp 

dist_data_DATA = parameters solver_test chkp_restart_test

