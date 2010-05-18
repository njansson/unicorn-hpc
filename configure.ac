 
AC_INIT([unicorn],[0.1.0-hpc])
AM_INIT_AUTOMAKE
AM_MAINTAINER_MODE
PKG_CHECK_MODULES([DOLFIN], [dolfin = 0.8.0-hpc])
PKG_CHECK_MODULES([PETSc], [petsc = 2.3.3])
CXXFLAGS="$DOLFIN_CFLAGS"
AC_PROG_CXX(`pkg-config --variable=compiler dolfin`)
AC_PROG_RANLIB
AC_CONFIG_FILES([Makefile src/Makefile	src/unicorn/Makefile unicorn.pc])
AC_OUTPUT
