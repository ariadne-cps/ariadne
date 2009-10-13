/***************************************************************************
 *            numeric_submodule.cc
 *
 *  Copyright 2008  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <iostream>
#include <iomanip>

#include "config.h"
#include "utilities.h"

#include "tribool.h"
#include "numeric.h"
#include "real.h"

#include <boost/python.hpp>

using namespace boost::python;

using namespace Ariadne;


namespace Ariadne {

bool definitely(bool b) { return b; }
bool possibly(bool b) { return b; }

void set_output_precision(uint p) { std::cout << std::setprecision(p); }


template<>
struct from_python_dict<Interval> {
    from_python_dict() { converter::registry::push_back(&convertible,&construct,type_id<Interval>()); }
    static void* convertible(PyObject* obj_ptr) {
        if (!PyDict_Check(obj_ptr) || len(extract<dict>(obj_ptr))!=1) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        boost::python::dict dct = boost::python::extract<boost::python::dict>(obj_ptr);
        boost::python::list lst=dct.items();
        assert(boost::python::len(lst)==1);
        void* storage = ((converter::rvalue_from_python_storage<Interval>*)data)->storage.bytes;
        new (storage) Interval(boost::python::extract<double>(lst[0][0]),boost::python::extract<double>(lst[0][1]));
        data->convertible = storage;
    }
};

template<>
struct from_python_list<Interval> {
    from_python_list() { converter::registry::push_back(&convertible,&construct,type_id<Interval>()); }
    static void* convertible(PyObject* obj_ptr) {
        if (!PyList_Check(obj_ptr) || len(extract<list>(obj_ptr))!=2) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        boost::python::list lst = boost::python::extract<boost::python::list>(obj_ptr);
        assert(boost::python::len(lst)==2);
        void* storage = ((converter::rvalue_from_python_storage<Interval>*)data)->storage.bytes;
        new (storage) Interval(boost::python::extract<double>(lst[0]),boost::python::extract<double>(lst[1]));
        data->convertible = storage;
    }
};

/*
struct interval_from_python_str {
    interval_from_python_str() { converter::registry::push_back(&convertible,&construct,type_id<Interval>()); }
    static void* convertible(PyObject* obj_ptr) { if (!PyString_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        std::string str = boost::python::extract<std::string>(obj_ptr);
        void* storage = ((converter::rvalue_from_python_storage<Interval>*)data)->storage.bytes;
        storage = new Interval(str);
        data->convertible = storage;
    }
};
*/
}


std::string __str__(tribool tb) {
  if(indeterminate(tb)) { return "Indeterminate"; }
  if(tb==true) { return "True"; }
  else if(tb==false) { return "False"; }
  else { return "Indeterminate"; }
}

std::string __repr__(tribool tb) {
  if(tb==true) { return "tribool(True)"; }
  else if(tb==false) { return "tribool(False)"; }
  else { return "tribool(Indeterminate)"; }
}

bool __nonzero__(tribool tb) { return tb==true; }
tribool indeterminate_const() { return indeterminate; }

void export_tribool() {

    class_<tribool> tribool_class("tribool",init<bool>());
    tribool_class.def(init<int>());
    tribool_class.def(init<tribool>());
    tribool_class.def("__eq__", (tribool(*)(tribool,tribool))(&boost::logic::operator==));
    tribool_class.def("__eq__", (tribool(*)(tribool,bool))(&boost::logic::operator==));
    tribool_class.def("__neq__", (tribool(*)(tribool,tribool))(&boost::logic::operator!=));
    tribool_class.def("__neq__", (tribool(*)(tribool,bool))(&boost::logic::operator!=));
    tribool_class.def("__and__", (tribool(*)(tribool,tribool))(&boost::logic::operator!=));
    tribool_class.def("__and__", (tribool(*)(tribool,bool))(&boost::logic::operator!=));
    tribool_class.def("__or__", (tribool(*)(tribool,tribool))(&boost::logic::operator!=));
    tribool_class.def("__or__", (tribool(*)(tribool,bool))(&boost::logic::operator!=));
    // WARNING: __not__ is not a special method!
    tribool_class.def("__not__", (tribool(*)(tribool))(&boost::logic::operator!));
    tribool_class.def("__nonzero__", (bool(*)(tribool))&__nonzero__);
    tribool_class.def("__str__", (std::string(*)(tribool))&__str__);
    tribool_class.def("__repr__", (std::string(*)(tribool))&__repr__);

    implicitly_convertible<bool,tribool>();

    def("indeterminate",(tribool(*)(void))&indeterminate_const);
    def("possibly",(bool(*)(bool))&possibly);
    def("possibly",(bool(*)(tribool))&possibly);
    def("definitely",(bool(*)(bool))&definitely);
    def("definitely",(bool(*)(tribool))&definitely);
    // no facility for wrapping C++ constants
    // def("Indeterminate",tribool_indeterminate_constant);

}

void export_float()
{
    def("set_output_precision", &set_output_precision);

    def("rec",(Float(*)(Float)) &Ariadne::rec);
    def("sqr",(Float(*)(Float)) &Ariadne::sqr);
    def("sqrt", (Float(*)(Float)) &Ariadne::sqrt);

    def("exp", (Float(*)(Float)) &Ariadne::exp);
    def("log", (Float(*)(Float)) &Ariadne::log);

    def("sin", (Float(*)(Float)) &Ariadne::sin);
    def("cos", (Float(*)(Float)) &Ariadne::cos);
    def("tan", (Float(*)(Float)) &Ariadne::tan);
    def("asin", (Float(*)(Float)) &Ariadne::asin);
    def("acos", (Float(*)(Float)) &Ariadne::acos);
    def("atan", (Float(*)(Float)) &Ariadne::atan);

}

#ifdef HAVE_GMPXX_H
void export_integer()
{
    class_<Integer> integer_class("Integer",init<int>());
    integer_class.def(boost::python::self_ns::str(self));
    //integer_class.def(boost::python::self_ns::repr(self));
    integer_class.def("__less__",(bool(*)(const mpz_class&, const mpz_class&)) &operator<);

    implicitly_convertible<int,Integer>();

}
#endif

#ifdef HAVE_GMPXX_H
void export_rational()
{
    class_<Rational> rational_class("Rational",init<int,int>());
    rational_class.def(init<int>());
    rational_class.def(boost::python::self_ns::str(self));
    //rational_class.def(boost::python::self_ns::repr(self));
    rational_class.def("__less__",(bool(*)(const Rational&, const Rational&)) &operator<);

    implicitly_convertible<int,Rational>();
    implicitly_convertible<double,Rational>();
    implicitly_convertible<Integer,Rational>();

}
#endif

void export_real()
{
    class_<Real> real_class("Real",init<Real>());
    real_class.def(init<Float>());
    real_class.def(init<Interval>());
    real_class.def(init<std::string>());
    real_class.def("radius", &Interval::radius);
    real_class.def(boost::python::self_ns::str(self));

    implicitly_convertible<int,Real>();
    implicitly_convertible<double,Real>();
    implicitly_convertible<Real,Interval>();
}


void export_interval()
{
    using boost::python::class_;
    using boost::python::init;
    using boost::python::self;
    using boost::python::return_value_policy;
    using boost::python::copy_const_reference;
    using boost::python::def;

    //typedef Interval (*IFUN)(const Interval&);
    //typedef Interval (*IZFUN)(const Interval&, int n);
    //typedef Interval (*INFUN)(const Interval&, unsigned int n);

    typedef Interval (*IFUN)(Interval);
    typedef Interval (*IIFUN)(Interval,Interval);
    typedef Interval (*IZFUN)(Interval, int n);
    typedef Interval (*INFUN)(Interval, unsigned int n);
    typedef bool (*IIPRED)(Interval,Interval);

    def("down",&down);
    def("up",&up);

    class_< Interval > interval_class("Interval");
    interval_class.def(init<double,double>());
    interval_class.def(init<double>());
    interval_class.def(init<Interval>());
#ifdef HAVE_GMPXX_H
    interval_class.def(init<Rational>());
    interval_class.def(init<Rational,Rational>());
#endif
    interval_class.def(+self);
    interval_class.def(-self);
    interval_class.def(self + self);
    interval_class.def(self - self);
    interval_class.def(self * self);
    interval_class.def(self / self);
    interval_class.def(self + double());
    interval_class.def(self - double());
    interval_class.def(self * double());
    interval_class.def(self / double());
    interval_class.def(double() + self);
    interval_class.def(double() - self);
    interval_class.def(double() * self);
    interval_class.def(double() / self);
    interval_class.def("lower", &Interval::lower, return_value_policy<copy_const_reference>());
    interval_class.def("upper", &Interval::upper, return_value_policy<copy_const_reference>());
    interval_class.def("midpoint", &Interval::midpoint);
    interval_class.def("radius", &Interval::radius);
    interval_class.def("width", &Interval::width);
    interval_class.def("contains", (bool(*)(Interval,Float)) &contains);
    interval_class.def(boost::python::self_ns::str(self));
    //interval_class.def(boost::python::self_ns::repr(self));

    implicitly_convertible<int,Interval>();
    implicitly_convertible<double,Interval>();
#ifdef HAVE_GMPXX_H
    implicitly_convertible<Rational,Interval>();
#endif

    from_python_dict<Interval>();
    from_python_list<Interval>();
    //from_python_str<Interval>();

    def("midpoint", &Interval::midpoint);
    def("radius", &Interval::radius);

    def("disjoint", (IIPRED) &disjoint);
    def("subset", (IIPRED) &subset);
    def("intersection", (IIFUN) &intersection);

    def("mag", (Float(*)(Interval)) &mag);

    def("med", (IFUN) &med);
    def("rad", (IFUN) &rad);
    def("diam", (IFUN) &diam);

    def("max", (IIFUN) &max);
    def("min", (IIFUN) &min);

    def("trunc", (INFUN) &trunc, "truncate to n binary digits");
    def("abs", (IFUN) &abs, "interval absolute value function");
    def("pow",  (IZFUN) &pow, "interval power function");
    def("sqr", (IFUN) &sqr, "interval square function");
    def("rec", (IFUN) &rec);
    def("sqrt", (IFUN) &sqrt);
    def("exp", (IFUN) &exp);
    def("log", (IFUN) &log);

    def("sin", (IFUN) &sin);
    def("cos", (IFUN) &cos);
    def("tan", (IFUN) &tan);
    def("asin", (IFUN) &asin);
    def("acos", (IFUN) &acos);
    def("atan", (IFUN) &atan);
}


void
numeric_submodule()
{
    export_tribool();
    export_float();
    export_interval();
    export_real();

#ifdef HAVE_GMPXX_H
    export_integer();
    export_rational();
#endif
}
