AC_PREREQ(2.57)
AC_INIT(ariadne, 0.1.0)
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

# Check for header files.
AC_STDC_HEADERS
AC_CHECK_HEADER([mpfr.h],,AC_MSG_ERROR([MPFR library is missing. You must install it.]))
AC_CHECK_HEADER([gmpxx.h],,AC_MSG_ERROR([GMP library is missing. You must install it.]))
AC_CHECK_HEADER([cblas.h],[AC_DEFINE([HAVE_CBLAS_H],[1],[Define to 1 if you have the <cblas.h> header file.])],AC_MSG_WARN([CBLAS library is missing. You cannot use Fortran BLAS routines.]))
AC_CHECK_HEADER([tblas/tblas.hpp],,AC_MSG_ERROR([TBLAS library is missing. You must install it.]))
AC_CHECK_HEADER([tlapack/tlapack.hpp],,AC_MSG_ERROR([TLAPACK library is missing. You must install it.]))
#AC_CHECK_HEADER([ppl.hh],[AC_DEFINE([HAVE_PPL_HH],[1],[Define to 1 if you have the <ppl.hh> header file.])],AC_MSG_WARN([PPL library is missing. Ariadne polyhedral routines cannot be checked against PPL routines.]))

# Check for Boost library
AC_CHECK_HEADER([boost/version.hpp],,AC_MSG_ERROR([Boost library is missing. You must install it.]))

# Find Python paths
AM_PATH_PYTHON(2.4,, :)
AC_MSG_NOTICE([PYTHON_VERSION=$PYTHON_PREFIX])
AC_MSG_NOTICE([PYTHON_VERSION=$PYTHON_VERSION])
AM_CHECK_PYTHON_HEADERS(,[AC_MSG_ERROR(could not find Python headers)])
CPPFLAGS="$CPPFLAGS $PYTHON_INCLUDES"

# Check for Boost python
#AC_CHECK_HEADER([boost/python.hpp],,AC_MSG_ERROR([Boost.Python library development files are missing. You must install them.]))

# Check for Doxygen
AC_PATH_PROG(DOXYGEN, doxygen)
if test -f "${DOXYGEN}"; then
    BUILD_DOC=yes;
else
    BUILD_DOC=no;
    AC_MSG_WARN([Doxygen is not installed. You can not generate Ariadne's documentation.])
fi


AC_MSG_CHECKING(for floating-point types to instantiate)
AC_ARG_ENABLE(float64,[AS_HELP_STRING([--enable-float64],[enable compilation of Float64 types (default is enabled)])],
              [],[enable_float64=yes])
AC_ARG_ENABLE(floatmp,[AS_HELP_STRING([--enable-floatmp],[enable compilation of FloatMP types (default is enabled)])],
              [],[enable_floatmp=yes])
AS_IF([test "x$enable_float64" != xno],[enable_float64="yes" AC_DEFINE([ENABLE_FLOAT64],[1],[Define to 1 if you want to enable Float64 compilation])],[])
AS_IF([test "x$enable_floatmp" != xno],[enable_floatmp="yes" AC_DEFINE([ENABLE_FLOATMP],[1],[Define to 1 if you want to enable FloatMP compilation])],[])
AC_MSG_RESULT([Float64: $enable_float64  FloatMP: $enable_floatmp])


AC_MSG_CHECKING([for float-point type in Python interface])
AC_ARG_ENABLE(python,[AS_HELP_STRING([--enable-python],[enable compilation of python floating-point types (default is Float64)])],
              [],[if test "x$enable_float64" != xno; then enable_python=Float64; elif test "x$enable_floatmp" != xno; then enable_python=FloatMP; else enable_python=no; fi])
#Make sure valid type is given for Python floating-point type
AS_IF([test "x$enable_python" != xno],
      [if test "x$enable_python" = "xfloat64" -o "x$enable_python" = "xFloat64"; then
         if test "x$enable_float64" = "xno"; then 
           AC_MSG_ERROR([Using Float64 for Python interface, but Float64 not instantiated]) 
         else
           enable_python=Float64; 
         fi
       elif test "x$enable_python" = "xfloatmp" -o "x$enable_python" = "xFloatMP"; then                                         
         if test "x$enable_floatmp" = "xno"; then 
           AC_MSG_ERROR([Using FloatMP for Python interface, but FloatMP not instantiated]) 
         else 
           enable_python=FloatMP; 
         fi
       else
         AC_MSG_ERROR([Unrecognised floating-point type for Python interface])
       fi
      ],[])      
AS_IF([test "x$enable_python" != xno],[AC_DEFINE_UNQUOTED([PYTHON_FLOAT],[$enable_python],[Define to the type used for floating-point numbers in Python])],[])
AC_MSG_RESULT([$enable_python])


#AC_CONFIG_FILES([Makefile include/Makefile src/Makefile src/base/Makefile src/numeric/Makefile src/linear_algebra/Makefile src/combinatoric/Makefile src/geometry/Makefile src/python/Makefile src/system/Makefile src/evaluation/Makefile test/Makefile])
AC_CONFIG_FILES([Makefile include/Makefile src/Makefile src/python/Makefile python/Makefile test/Makefile])
AC_OUTPUT
