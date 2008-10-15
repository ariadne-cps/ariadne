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
 
#include "numeric.h"

#include <boost/python.hpp>

using namespace boost::python;

using namespace Ariadne;

void read(Interval& ivl, const boost::python::object& obj) {
  boost::python::extract<boost::python::list> lx(obj);
  boost::python::extract<boost::python::tuple> tx(obj);
  boost::python::extract<double> dx(obj);
  boost::python::extract<Interval> ix(obj);

  if(lx.check()) {
    boost::python::list lst=lx();
    assert(boost::python::len(lst)==2);
    double l=boost::python::extract<double>(lst[0]);
    double u=boost::python::extract<double>(lst[1]);
    ivl=Interval(l,u);
  } else if(tx.check()) {
    boost::python::tuple tpl=tx();
    assert(boost::python::len(tpl)==2);
    double l=boost::python::extract<double>(tpl[0]);
    double u=boost::python::extract<double>(tpl[1]);
    ivl=Interval(l,u);
  } else if(dx.check()) {
    ivl=static_cast<Interval>(dx());
  } else if(ix.check()) {
    ivl=ix();
  }
}



void export_float() 
{
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
    interval_class.def(boost::python::self_ns::str(self));

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
  export_float();
  export_interval();
}
