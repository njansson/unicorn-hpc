AM_LDFLAGS = @DOLFIN_LIBS@
AM_CPPFLAGS = -I@top_builddir@/src
check_PROGRAMS = minimal
minimal_SOURCES = minimal.cpp
LDADD = $(top_builddir)/src/libunicorn.la
TESTS = minimal_solver_test

dist_data_DATA = parameters minimal_solver_test usquare.xml