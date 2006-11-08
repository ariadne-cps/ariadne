AC_PREREQ(2.57)
AC_INIT(ariadne, 0.0.0)
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src])
AC_CONFIG_HEADER([config.h])

dnl Checks for programs.
AC_PROG_CXX

dnl Checks for libraries.
AC_LANG(C++)

# Disable static libraries
AC_DISABLE_STATIC

# Checks for Libtool
AC_PROG_LIBTOOL

# Check type sizes
AC_CHECK_SIZEOF([int])
AC_CHECK_SIZEOF([size_t])
AC_CHECK_SIZEOF([long int])

dnl Checks for header files.
AC_STDC_HEADERS
AC_CHECK_HEADER([mpfr.h],,AC_MSG_ERROR([MPFR library is missing. You must install it.]))
AC_CHECK_HEADER([gmpxx.h],,AC_MSG_ERROR([GMP library is missing. You must install it.]))
AC_CHECK_HEADER([tblas/tblas.hpp],,AC_MSG_ERROR([TBLAS library is missing. You must install it.]))
AC_CHECK_HEADER([tlapack/tlapack.hpp],,AC_MSG_ERROR([TLAPACK library is missing. You must install it.]))
AC_CHECK_HEADER([ppl.hh],,AC_MSG_ERROR([PPL library is missing. You must install it.]))

CPPFLAGS="$CPPFLAGS -I/usr/include/python2.4"
AC_CHECK_HEADER([boost/version.hpp],,AC_MSG_ERROR([Boost library is missing. You must install it.]))
AC_CHECK_HEADER([boost/python.hpp],,AC_MSG_ERROR([Boost.Python library development files are missing. You must install them.]))

#AC_MSG_CHECKING(for Doxygen)
AC_PATH_PROG(DOXYGEN, doxygen)
if test -f "${DOXYGEN}"; then
    BUILD_DOC=yes;
else
    BUILD_DOC=no;
    AC_MSG_WARN([Doxygen is not installed. You can not generate Ariadne's documentation.])
fi

# Find Python paths
AM_PATH_PYTHON(2.4,, :)

#AC_CONFIG_FILES([Makefile include/Makefile src/Makefile src/base/Makefile src/numeric/Makefile src/linear_algebra/Makefile src/combinatoric/Makefile src/geometry/Makefile src/python/Makefile src/system/Makefile src/evaluation/Makefile test/Makefile])
AC_CONFIG_FILES([Makefile include/Makefile src/Makefile src/python/Makefile python/Makefile test/Makefile])
AC_OUTPUT
