# Declare versions
AC_PREREQ(2.57)
AC_INIT(ariadne, 0.4)
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CXX

# Checks for libraries.
AC_LANG(C++)

# Disable static libraries
AC_DISABLE_STATIC

# Checks for Libtool
AC_PROG_LIBTOOL

# Check type sizes
AC_CHECK_SIZEOF([int])
AC_CHECK_SIZEOF([size_t])
AC_CHECK_SIZEOF([long int])

# Check for standard header files.
AC_STDC_HEADERS

# Check for Boost library (required)
AC_CHECK_HEADER([boost/version.hpp],,AC_MSG_ERROR([Boost library is missing. You must install it.]))

# Check for GMP (optional)
AC_CHECK_HEADER([gmpxx.h],[AC_DEFINE([HAVE_GMPXX_H],[1],[Define to 1 if you have the <gmpxx.h> header file.])],AC_MSG_WARN([GNU Multiple Precision (GMP) library development files are missing. You cannot use Integer or Rational numbers.]))

# Check for MPFR (optional)
AC_CHECK_HEADER([mpfr.h],[AC_DEFINE([HAVE_MPFR_H],[1],[Define to 1 if you have the <mpfr.h> header file.])],AC_MSG_WARN([MPFR library development files are missing. You cannot use multiple-precision floating-point numbers.]))

# Check for GTK (optional)
PKG_CHECK_MODULES(GTK, gtk+-2.0 >= 2.2, [AC_DEFINE([HAVE_CAIRO_H],[1],[Define to 1 if you have the <gtk.h> header file.])],AC_MSG_WARN([GTK library development files are missing. You cannot use graphics windows.]) )

# Check for Cairo (optional)
PKG_CHECK_MODULES(CAIRO, cairo, [AC_DEFINE([HAVE_CAIRO_H],[1],[Define to 1 if you have the <cairo.h> header file.])],AC_MSG_WARN([Cairo library development files are missing. You cannot use graphics facilities.]))

# Find Python paths
AM_PATH_PYTHON(2.5,, :)
AC_MSG_NOTICE([PYTHON_VERSION=$PYTHON_PREFIX])
AC_MSG_NOTICE([PYTHON_VERSION=$PYTHON_VERSION])
AM_CHECK_PYTHON_HEADERS(,[AC_MSG_ERROR(could not find Python headers)])
CPPFLAGS="$CPPFLAGS $PYTHON_INCLUDES"

# Check for Boost python library (optional)
AC_CHECK_HEADER([boost/python.hpp],[AC_DEFINE([HAVE_BOOST_PYTHON_HPP],[1],[Define to 1 if you have the <boost/python.hpp> header file.])],AC_MSG_WARN([Boost Python library development files are missing. You cannot compile the Python interface.]))

# Check for Doxygen (optional)
AC_PATH_PROG(DOXYGEN, doxygen)
if test -f "${DOXYGEN}"; then
    BUILD_DOC=yes;
else
    BUILD_DOC=no;
    AC_MSG_WARN([Doxygen is not installed. You can not generate Ariadne's documentation.])
fi


# Check whether to use multiple-precision floating-point (default no)
AC_MSG_CHECKING(whether to use multiple-precision floating-point numbers)
AC_ARG_ENABLE(mp,[AS_HELP_STRING([--enable-mp],[enable compilation of FloatMP types (default is disabled)])],
              [AC_DEFINE([ENABLE_MULTIPLE_PRECISION],[1],[Define to 1 if you want to use multiple-precision floating-point numbers.])],[AC_MSG_RESULT([no])],[AC_MSG_RESULT([yes])])



# Declare which Makefiles need to be built
AC_CONFIG_FILES([Makefile include/Makefile src/Makefile wrap/Makefile python/Makefile test/Makefile examples/Makefile])

# Output results of autoconf tests
AC_OUTPUT


