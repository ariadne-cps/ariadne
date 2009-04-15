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
 
#include "config.h"
#include "utilities.h"

#include "tribool.h"
#include "numeric.h"

#include <boost/python.hpp>

using namespace boost::python;

using namespace Ariadne;


namespace Ariadne {
bool definitely(bool b) { return b; }
bool possibly(bool b) { return b; }
}

const char * tribool_c_str(tribool tb) {
  if(indeterminate(tb)) { return "Indeterminate"; }
  if(tb==true) { return "True"; }
  else if(tb==false) { return "False"; }
  else { return "Indeterminate"; }
}

const char * tribool_c_repr(tribool tb) {
  if(tb==true) { return "tribool(True)"; }
  else if(tb==false) { return "tribool(False)"; }
  else { return "tribool(Indeterminate)"; }
}

tribool tribool_not(tribool tb) { return !tb; }
bool tribool_nonzero(tribool tb) { return tb==true; }
tribool tribool_indeterminate() { return indeterminate; }

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
    tribool_class.def("__nonzero__", &tribool_nonzero);
    tribool_class.def("__str__", &tribool_c_str);
    tribool_class.def("__repr__", &tribool_c_repr);
  
    def("indeterminate",(tribool(*)(void))&tribool_indeterminate);
    def("possibly",(bool(*)(bool))&possibly);
    def("possibly",(bool(*)(tribool))&possibly);
    def("definitely",(bool(*)(bool))&definitely);
    def("definitely",(bool(*)(tribool))&definitely);
    // no facility for wrapping C++ constants
    // def("Indeterminate",tribool_indeterminate_constant);

}

void export_float() 
{
}

#ifdef HAVE_GMPXX_H
void export_rational()
{
    class_<Rational> rational_class("Rational",init<int,int>());
    rational_class.def(boost::python::self_ns::str(self));
    rational_class.def("__less__",(bool(*)(const Rational&, const Rational&)) &operator<);
}
#endif

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
    typedef Interval (*IZFUN)(Interval, int n);
    typedef Interval (*INFUN)(Interval, unsigned int n);
    
    def("down",&down);
    def("up",&up);

    class_< Interval > interval_class("Interval");
    interval_class.def(init<double,double>());
    interval_class.def(init<double>());
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
    interval_class.def("__repr__",&__repr__<Interval>);
    interval_class.def(boost::python::self_ns::str(self));

    def("midpoint", &Interval::midpoint);
    def("radius", &Interval::radius);
    def("disjoint", &disjoint);
    def("subset", &subset);
    def("intersection", &intersection);

    def("mag", (Float(*)(Interval)) &mag);

    def("med", (IFUN) &med);
    def("rad", (IFUN) &rad);
    def("diam", (IFUN) &diam);

    def("trunc", (INFUN) &trunc, "truncate to n binary digits");
    def("abs", (IFUN) &abs, "interval absolute value function");
    def("pow",  (IZFUN) &pow, "interval power function");
    def("sqr", (IFUN) &sqr, "interval square function");
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

#ifdef HAVE_GMPXX_H
    export_rational();
#endif
}
