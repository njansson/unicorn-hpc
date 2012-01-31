 
AC_INIT([unicorn],[0.1.3-hpc])
AC_DEFINE_UNQUOTED(UNICORN_VERSION,["${PACKAGE_VERSION}"],[Version of Unicorn])
AM_INIT_AUTOMAKE
AM_MAINTAINER_MODE
AC_LANG(C++)
AC_CONFIG_MACRO_DIR([m4])
LT_INIT([disable-shared])
PKG_CHECK_MODULES([DOLFIN], [dolfin = 0.8.2-hpc])
CXXFLAGS="$CXXFLAGS $DOLFIN_CFLAGS"
AC_PROG_CXX(`pkg-config --variable=compiler dolfin`)
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

AC_HEADER_STDBOOL
AC_CHECK_HEADERS([stdlib.h sys/time.h unistd.h])
AC_CHECK_FUNCS([gettimeofday memset mkdir],,
			     AC_MSG_ERROR([Vital functions are missing]))
AC_FUNC_ERROR_AT_LINE


# Try to find a suitable BLAS
AC_CHECK_HEADERS([cblas.h gsl_cblas.h sunperf.h scsl_cblas.h],,)
have_blas=no
AC_SEARCH_LIBS([sgemm], [sunperf scs blas sci goto2],have_blas=yes, have_blas=no)
if test $have_blas = no; then
   AC_SEARCH_LIBS([cblas_sgemm], [sunperf scs blas gslcblas sci goto2],have_blas=yes, have_blas=no)
   if test $have_blas = no; then
      AC_CHECK_LIB(cblas, cblas_sgemm, [have_blas=yes, LIBS="-lcblas -latlas -lm $LIBS"],[have_blas=no],[-latlas -lm])
   fi
   
   # Fallback to FORTRAN BLAS
   if test $have_blas = no; then
      AC_SEARCH_LIBS([sgemm_], [sunperf scs_mp scs blas sci goto2_itanium2p-r1.13],have_blas=yes, have_blas=no)
      if test $have_blas = yes; then
         AC_DEFINE(HAVE_F77_BLAS, [1], [Have FORTRAN BLAS])
      fi
   fi

fi

if test $have_blas = no; then
   AC_MSG_ERROR([suitable BLAS not found])
fi


AC_CONFIG_HEADERS([src/unicorn/unicorn_config.h])
AC_CONFIG_FILES([Makefile src/Makefile test/Makefile \
			  test/unicorn/Makefile \
			  test/icns/Makefile unicorn.pc])
AC_OUTPUT

# Display help text
echo \
"
---------------------------------------------------------

Configuration of Unicorn ${PACKAGE_VERSION} finished.

Prefix: '${prefix}'.
Compiler: ${CXX} ${CXXFLAGS} ${CPPFLAGS}
Libs: ${LD} ${LDFLAGS} ${LIBS}

---------------------------------------------------------
"
