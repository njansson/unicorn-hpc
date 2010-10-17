 
AC_INIT([unicorn],[0.1.1-hpc])
AM_INIT_AUTOMAKE
AM_MAINTAINER_MODE
AC_LANG(C++)
AC_CONFIG_MACRO_DIR([m4])
LT_INIT([disable-shared])
PKG_CHECK_MODULES([DOLFIN], [dolfin = 0.8.1-hpc])
CXXFLAGS="$DOLFIN_CFLAGS"
AC_PROG_CXX(`pkg-config --variable=compiler dolfin`)
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_C_INLINE

AC_HEADER_STDBOOL
AC_CHECK_HEADERS([stdlib.h sys/time.h unistd.h])
AC_CHECK_FUNCS([gettimeofday memset mkdir],,
			     AC_MSG_ERROR([Vital functions are missing]))
AC_FUNC_ERROR_AT_LINE

AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_FILES([Makefile src/Makefile unicorn.pc])
AC_OUTPUT

