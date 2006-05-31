AC_PREREQ(2.57)
AC_INIT(ariadne, 0.0.0)
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src])
AC_CONFIG_HEADER([config.h])

dnl Checks for programs.
AC_PROG_CXX

dnl Checks for libraries.
AC_LANG(C++)

dnl Checks for header files.
AC_STDC_HEADERS
AC_CHECK_HEADER([gmpxx.h],,AC_MSG_ERROR([GMP library is missing. You must install it.]))

CPPFLAGS="$CPPFLAGS -I/usr/include/python2.4"
AC_CHECK_HEADER([boost/version.hpp],,AC_MSG_ERROR([BOOST library is missing. You must install it.]))
AC_CHECK_HEADER([boost/python.hpp],,AC_MSG_ERROR([Boost.Python library development files are missing. You must install them.]))

#AC_MSG_CHECKING(for Doxygen)
AC_PATH_PROG(DOXYGEN, doxygen)
if test -f "${DOXYGEN}"; then
    BUILD_DOC=yes;
else
    BUILD_DOC=no;
    AC_MSG_WARN([Doxygen is not installed. You can not generate Ariadne's documentation.])
fi


AC_CONFIG_FILES([Makefile src/Makefile src/base/Makefile src/evaluation/Makefile src/geometry/Makefile src/linear_algebra/Makefile src/numeric/Makefile src/python/Makefile])
AC_OUTPUT
