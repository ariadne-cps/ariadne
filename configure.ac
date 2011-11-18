# Declare versions
AC_PREREQ(2.57)
AC_INIT(ariadne, 0.4)
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])

# Save the user-given value of CXXFLAGS for later use
CXXFLAGS_EARLY=${CXXFLAGS}

# Checks for programs.
AC_PROG_CXX

# this macro is used to get the arguments supplied
# to the configure script (./configure --enable-debug)
# Check if we have enable debug support.
AC_MSG_CHECKING(whether to enable debugging)
AC_ARG_ENABLE(debug, [AS_HELP_STRING([--enable-debug=(yes/no)], [turn on debugging (default no)])],
              enable_debug=${enableval}, enable_debug="no")
# If enable_debug is yes, then put the correct CXXFLAGS
if test "x$enable_debug" = "xyes"; then
    # If user provided a custom value for CXXFLAGS, do no change it
    if test "x$CXXFLAGS_EARLY" = "x"; then
        CXXFLAGS="-pg -g3 -O0"
    fi
    AC_MSG_RESULT(yes)
else
    AC_MSG_RESULT(no)
fi

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

# Check for standard fenv.h (optional)
AC_CHECK_HEADER([fenv.h],[AC_DEFINE([HAVE_FENV_H],[1],[Define to 1 if you have the <fenv.h> header file.])],[AC_MSG_WARN([Standard C header file fenv.h is missing. Only assembly-language rounding is supported.])])

# Check for Boost library (required)
BOOST_REQUIRE([1.36.0])
#AC_CHECK_HEADER([boost/version.hpp],,[AC_MSG_ERROR([Boost library development files are missing. You must install them.])])

# Check for GMP (optional)
AC_CHECK_HEADER([gmpxx.h],[AC_DEFINE([HAVE_GMPXX_H],[1],[Define to 1 if you have the <gmpxx.h> header file.])],[AC_MSG_WARN([GNU Multiple Precision (GMP) library development files are missing. You cannot use Integer or Rational numbers.])])

# Check for MPFR (optional)
AC_CHECK_HEADER([mpfr.h],[AC_DEFINE([HAVE_MPFR_H],[1],[Define to 1 if you have the <mpfr.h> header file.])],[AC_MSG_WARN([MPFR library development files are missing. You cannot use multiple-precision floating-point numbers.])])

# Check for GTK (optional)
PKG_CHECK_MODULES(GTK, gtk+-2.0 >= 2.2, [AC_DEFINE([HAVE_GTK_H],[1],[Define to 1 if you have the <gtk.h> header file.])],[AC_MSG_WARN([GTK+-2.0 library development files are missing. You cannot use graphics windows.])])

# Check for Cairo (optional)
PKG_CHECK_MODULES(CAIRO, cairo, [AC_DEFINE([HAVE_CAIRO_H],[1],[Define to 1 if you have the <cairo.h> header file.])],[AC_MSG_WARN([Cairo library development files are missing. You cannot use graphics facilities.])])

# this macro is used to get the arguments supplied
# to the configure script (./configure --enable-python)
# Check if we have enable python support.
AC_MSG_CHECKING(whether to enable python)
AC_ARG_ENABLE(python, [AS_HELP_STRING([--enable-python=(yes/no)], [turn on the python library (default yes)])],
              enable_python=${enableval}, enable_python="yes")
# If enable_python is yes, check for the python libraries
if test "x$enable_python" = "xyes"; then
    AC_MSG_RESULT(yes)
    # Find Python paths (optional)
    AM_PATH_PYTHON(2.5,, :)
    AC_MSG_NOTICE([PYTHON_VERSION=$PYTHON_VERSION])
    AM_CHECK_PYTHON_HEADERS([PYTHON_LIB=-lpython$PYTHON_VERSION ],[AC_MSG_WARN([could not find Python headers. You cannot compile the Python interface.])])
    AC_SUBST(PYTHON_LIB)
    # Add Python include directory to include path since otherwise check for Boost Python fails
    CPPFLAGS+=" $PYTHON_INCLUDES "
    # Check for Boost python library (optional)
    AC_CHECK_HEADER([boost/python.hpp],[AC_DEFINE([HAVE_BOOST_PYTHON_HPP],[1],[Define to 1 if you have the <boost/python.hpp> header file.]) HAVE_BOOST_PYTHON_HPP=true],[AC_MSG_WARN([Boost Python library development files are missing. You cannot compile the Python interface.])])
    if test x$HAVE_BOOST_PYTHON_HPP = xtrue ; then
        BOOST_PYTHON
    fi
else
    AC_MSG_RESULT(no)
fi

# Check for GMP library (optional)
AC_CHECK_LIB([gmp],[main])

# Check for GMPXX library (optional)
AC_CHECK_LIB([gmpxx],[main])

# Check for MPFR library (optional)
AC_CHECK_LIB([mpfr],[main])

# Check for Boost Serialization library (optional)
# this macro is used to get the arguments supplied
# to the configure script (./configure --enable-serialization)
# Check if we have enable serialization support.
AC_MSG_CHECKING(whether to enable serialization)
AC_ARG_ENABLE(serialization, [AS_HELP_STRING([--enable-serialization=(yes/no)], [turn on serialization (default no)])],
              enable_serialization=${enableval}, enable_serialization="no")
# If enable_serialization is yes, then check for Boost Serialization library
if test "x$enable_serialization" = "xyes"; then
    AC_MSG_RESULT(yes)
    BOOST_SERIALIZATION
    AC_DEFINE([ARIADNE_ENABLE_SERIALIZATION],[1],[Define to 1 if you want to use the Boost serialization library.])
else
    AC_MSG_RESULT(no)
fi 

AM_CONDITIONAL([WANT_PYTHON_INTERFACE], [test x$HAVE_BOOST_PYTHON_HPP = xtrue])
AM_CONDITIONAL([WANT_SERIALIZATION], [test "x$enable_serialization" = "xyes"])

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
              [AC_DEFINE([ENABLE_MULTIPLE_PRECISION],[1],[Define to 1 if you want to use multiple-precision floating-point numbers.])],
              [AC_MSG_RESULT([no])],[AC_MSG_RESULT([yes])])

# this macro is used to get the arguments supplied
# to the configure script (./configure --with-static-exec)
# Check if we have to link executables statically
# (needed to allow debugging in Eclipse).
AC_MSG_CHECKING(whether to to link executables statically)
AC_ARG_ENABLE(static-exec, [AS_HELP_STRING([--enable-static-exec(=yes/no)],[link statically executable files (default no)])],
            static_exec=${enableval}, static_exec="no")
# If static_exec is yes, then put the correct EXEC_LDFLAGS
if test "x$static_exec" = "xyes"; then
    EXEC_LDFLAGS="-static"
    AC_MSG_RESULT(yes)
else
    EXEC_LDFLAGS=""
    AC_MSG_RESULT(no)
fi

# C++ flags to be used for all files
#   Add -fno-strict-aliasing because of code in Expansion class which uses
#   reinterpret_cast in a way that breaks strict anti-aliasing rules.
#   Add -frounding-math to allow correctly-rounded arithmetic.
ARIADNE_CXXFLAGS="${ARIADNE_CXXFLAGS} -Wall -fno-strict-aliasing -frounding-math ${BOOST_CPPFLAGS}"

AC_MSG_NOTICE([LIBS=$LIBS])
AC_MSG_NOTICE([GTK_LIBS=$GTK_LIBS])
AC_MSG_NOTICE([CAIRO_LIBS=$CAIRO_LIBS])
AC_MSG_NOTICE([PYTHON_LIB=$PYTHON_LIB])
AC_MSG_NOTICE([BOOST_PYTHON_LIB=$BOOST_PYTHON_LIB])
AC_MSG_NOTICE([HAVE_BOOST_PYTHON_HPP='$HAVE_BOOST_PYTHON_HPP'])
AC_MSG_NOTICE([HAVE_LIBBOOST_PYTHON='$HAVE_LIBBOOST_PYTHON'])
AC_MSG_NOTICE([WANT_PYTHON_INTERFACE=$WANT_PYTHON_INTERFACE])
AC_MSG_NOTICE([ARIADNE_CXXFLAGS=$ARIADNE_CXXFLAGS])
AC_MSG_NOTICE([EXEC_LDFLAGS=$EXEC_LDFLAGS])

AC_SUBST(BOOST_PYTHON_LIB)
AC_SUBST(EXEC_LDFLAGS)
AC_SUBST(ARIADNE_CXXFLAGS)

# Declare which Makefiles need to be built
AC_CONFIG_FILES([Makefile include/Makefile src/Makefile wrap/Makefile python/Makefile test/Makefile examples/Makefile])

# Output results of autoconf tests
AC_OUTPUT


