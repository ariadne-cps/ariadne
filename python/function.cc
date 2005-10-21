// g++ -fpic -I../include/ -I/usr/include/python2.4 -c interval.cc
// g++ -shared interval.o -o interval.so -lboost_python -lmpfr -lgmp  -lsuperlu -lblas

#include <iostream>
#include <boost/python.hpp>
#include <gmpxx.h>

#include "function.h"

typedef Ariadne::Dyadic Real;

BOOST_PYTHON_MODULE(function)
{
    using boost::python::class_;
    using boost::python::init;
    using boost::python::self;
    using boost::python::return_value_policy;
    using boost::python::copy_const_reference;
    using boost::python::def;

    def("approx_exp", Ariadne::approx_exp<Real>, "approximate exponential function (maximum error e)" );
    def("approx_cos", Ariadne::approx_exp<Real>, "approximate sine function (maximum error e)" );
    def("approx_sin", Ariadne::approx_exp<Real>, "approximate cosine function (maximum error e)" );
}
