/***************************************************************************
 *            calculus_submodule.cc
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
 
#include "taylor_variable.h"

#include <boost/python.hpp>
using namespace boost::python;

using namespace Ariadne;

#include "utilities.h"

namespace Ariadne {

    
void read(MultiIndex& j, const boost::python::object& obj) {
    array<int> ary;
    read_tuple_array(ary,obj);
    j=MultiIndex(ary.size(),ary.begin());
}


void read(Polynomial<Float>& p, const boost::python::object& obj) {
    Float c;
    MultiIndex j(0);
    boost::python::dict dct=extract<boost::python::dict>(obj);
    boost::python::list lst=dct.items();
    for(int i=0; i!=len(lst); ++i) {
        boost::python::tuple tup=extract<boost::python::tuple>(lst[i]);
        read(j,tup[0]);
        read(c,tup[1]);
        if(i==0) { 
            p=Polynomial<Float>(j.size());
        }
        p[j]=c;
    }
}



void read(TaylorVariable& tv, const boost::python::object& obj1, const boost::python::object& obj2) {
    Polynomial<Float> p;
    Float e;
    read(p,obj1);
    read(e,obj2); 
    ARIADNE_ASSERT(e>=0);
    tv=TaylorVariable(p,e);
}
   

void read(TaylorVariable& tv, const boost::python::object& obj1, const boost::python::object& obj2,
          const boost::python::object& obj3, const boost::python::object& obj4) 
{
    uint as;
    uint deg;
    array<Float> dat;
    Float err;
    read(as,obj1);
    read(deg,obj2);
    read_list_array(dat,obj3);
    read(err,obj4);
    const Float* ptr=dat.begin();
    ARIADNE_ASSERT(err>=0);
    tv=TaylorVariable(as,deg,ptr,err);
}

void read(TaylorVariable& tv, const boost::python::object& obj);

void read(TaylorVariable& tv, const boost::python::tuple& tup) {
    if(len(tup)==1) {
        read(tv,tup[0]);
    } else if(len(tup)==2) {
        read(tv,tup[0],tup[1]);
    } else if(len(tup)==4) {
        read(tv,tup[0],tup[1],tup[2],tup[3]);
    }
}


void read(TaylorVariable& tv, const boost::python::object& obj) {
    if(check(extract<boost::python::tuple>(obj))) {
         read(tv,extract<boost::python::tuple>(obj)); 
    } else {
        Polynomial<Float> p;
        Float e=0;
        read(p,obj);
        tv=TaylorVariable(p,e);
    } 
}
   
TaylorVariable*
make_taylor_variable(const uint& as, const uint& deg, const boost::python::object& aobj, const boost::python::object& eobj) 
{
    array<Float> dat;
    Float err;
    read_list_array(dat,aobj);
    const Float* ptr=dat.begin();
    read(err,eobj);
    ARIADNE_ASSERT(err>=0);
    return new TaylorVariable(as,deg,ptr,err);
}

Vector<TaylorVariable>
make_taylor_variables(const Vector<Interval>& x) 
{
    Vector<TaylorVariable> result(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        result[i]=TaylorVariable::variable(x.size(),midpoint(x[i]),i);
    }
    return result;
}

   


boost::python::tuple
python_split(const TaylorVariable& x, uint j)
{
    std::pair<TaylorVariable,TaylorVariable> res=split(x,j);
    return boost::python::make_tuple(res.first,res.second);
}



template<> std::string __repr__(const TaylorVariable& tv) {
    std::stringstream ss;
    ss << "TaylorVariable(" << tv.expansion() << "," << tv.error() << ")";
    return ss.str();
} 

TaylorVariable neg(const TaylorVariable&);
TaylorVariable rec(const TaylorVariable&);
TaylorVariable sqr(const TaylorVariable&);
TaylorVariable pow(const TaylorVariable&, int);
TaylorVariable sqrt(const TaylorVariable&);
TaylorVariable exp(const TaylorVariable&);
TaylorVariable log(const TaylorVariable&);
TaylorVariable sin(const TaylorVariable&);
TaylorVariable cos(const TaylorVariable&);
TaylorVariable tan(const TaylorVariable&);

} // namespace Ariadne

void export_taylor_variable() 
{
    typedef uint N;
    typedef double D;
    typedef Float R;
    typedef Interval I;
    typedef MultiIndex A;
    typedef Vector<Float> RV;
    typedef Vector<Interval> IV;
    typedef Matrix<Float> RMx;
    typedef TaylorVariable T;
    typedef Vector<TaylorVariable> TV;

    class_<T> taylor_variable_class("TaylorVariable");
    taylor_variable_class.def("__init__", make_constructor(&make<TaylorVariable>) );
    taylor_variable_class.def("__init__", make_constructor(&make2<TaylorVariable>) );
    taylor_variable_class.def("__init__", make_constructor(&make_taylor_variable) );
    taylor_variable_class.def( init< uint >());
    taylor_variable_class.def( init< TaylorVariable >());
    taylor_variable_class.def("error", (const R&(T::*)()const) &T::error, return_value_policy<copy_const_reference>());
    taylor_variable_class.def("domain", &TaylorVariable::domain);
    taylor_variable_class.def("range", &TaylorVariable::range);
    taylor_variable_class.def("__getitem__", &get_item<T,A,R>);
    taylor_variable_class.def("__setitem__",&set_item<T,A,D>);
    taylor_variable_class.def(-self);
    taylor_variable_class.def(self+self);
    taylor_variable_class.def(self-self);
    taylor_variable_class.def(self*self);
    taylor_variable_class.def(self/self);
    taylor_variable_class.def(self+R());
    taylor_variable_class.def(self-R());
    taylor_variable_class.def(self*R());
    taylor_variable_class.def(self/R());
    taylor_variable_class.def(R()+self);
    taylor_variable_class.def(R()-self);
    taylor_variable_class.def(R()*self);
    taylor_variable_class.def(R()/self);
    taylor_variable_class.def(self+=R());
    taylor_variable_class.def(self-=R());
    taylor_variable_class.def(self*=R());
    taylor_variable_class.def(self/=R());
    taylor_variable_class.def(self+=self);
    taylor_variable_class.def(self-=self);
    taylor_variable_class.def(self_ns::str(self));
    taylor_variable_class.def("truncate", (TaylorVariable&(TaylorVariable::*)(uint)) &TaylorVariable::truncate,return_value_policy<reference_existing_object>());
    taylor_variable_class.def("sweep", (TaylorVariable&(TaylorVariable::*)(double))&TaylorVariable::sweep,return_value_policy<reference_existing_object>());
    taylor_variable_class.def("truncate", (TaylorVariable&(TaylorVariable::*)()) &TaylorVariable::truncate,return_value_policy<reference_existing_object>());
    taylor_variable_class.def("sweep", (TaylorVariable&(TaylorVariable::*)())&TaylorVariable::sweep,return_value_policy<reference_existing_object>());
    taylor_variable_class.def("clean", (TaylorVariable&(TaylorVariable::*)()) &TaylorVariable::clean,return_value_policy<reference_existing_object>());
    taylor_variable_class.def("set_maximum_degree",&TaylorVariable::set_maximum_degree);
    taylor_variable_class.def("set_sweep_threshold",&TaylorVariable::set_sweep_threshold);
    taylor_variable_class.def("maximum_degree",&TaylorVariable::maximum_degree);
    taylor_variable_class.def("sweep_threshold",&TaylorVariable::sweep_threshold);
    taylor_variable_class.def("set_default_maximum_degree",&TaylorVariable::set_default_maximum_degree);
    taylor_variable_class.def("set_default_sweep_threshold",&TaylorVariable::set_default_sweep_threshold);
    taylor_variable_class.def("default_maximum_degree",&TaylorVariable::default_maximum_degree);
    taylor_variable_class.def("default_sweep_threshold",&TaylorVariable::default_sweep_threshold);
    taylor_variable_class.def("evaluate", (Interval(TaylorVariable::*)(const Vector<Float>&)const) &TaylorVariable::evaluate);
    taylor_variable_class.def("evaluate", (Interval(TaylorVariable::*)(const Vector<Interval>&)const) &TaylorVariable::evaluate);
    taylor_variable_class.staticmethod("set_default_maximum_degree");
    taylor_variable_class.staticmethod("set_default_sweep_threshold");
    taylor_variable_class.staticmethod("default_maximum_degree");
    taylor_variable_class.staticmethod("default_sweep_threshold");
    def("split", &python_split);


    taylor_variable_class.def("__repr__",&__repr__<T>);

    taylor_variable_class.def("constant",(T(*)(N, const R&))&T::constant);
    taylor_variable_class.def("variable",(T(*)(N, const R&, uint))&T::variable);
    taylor_variable_class.def("variables",&make_taylor_variables);
  
    taylor_variable_class.staticmethod("constant");
    taylor_variable_class.staticmethod("variable");
    taylor_variable_class.staticmethod("variables");

    def("join", (TV(*)(const TV&,const TV&)) &join);
    def("join", (TV(*)(const TV&,const T&)) &join);

    def("unscale", (T(*)(const T&,const I&)) &unscale);
    def("unscale", (TV(*)(const TV&,const IV&)) &unscale);
 
    def("scale", (T(*)(const T&,const I&)) &scale);
    def("scale", (TV(*)(const TV&,const IV&)) &scale);

    def("evaluate",(IV(*)(const TV&,const IV&)) &evaluate);
    def("evaluate",(I(*)(const T&,const IV&)) &evaluate);

    def("compose",(TV(*)(const TV&,const IV&,const TV&)) &compose);
    def("compose",(T(*)(const T&,const IV&,const TV&)) &compose);
    def("compose",(T(*)(const T&,const I&,const T&)) &compose);

    def("max",(T(*)(const T&,const T&))&max);
    def("min",(T(*)(const T&,const T&))&min);
    def("abs",(T(*)(const T&))&abs);
    def("neg",(T(*)(const T&))&neg);
    def("rec",(T(*)(const T&))&rec);
    def("sqr",(T(*)(const T&))&sqr);
    def("pow",(T(*)(const T&, int))&pow);

    def("sqrt", (T(*)(const T&))&sqrt);
    def("exp", (T(*)(const T&))&exp);
    def("log", (T(*)(const T&))&log);
    def("sin", (T(*)(const T&))&sin);
    def("cos", (T(*)(const T&))&cos);
    def("tan", (T(*)(const T&))&tan);

    class_< Vector<TaylorVariable> >("TaylorVariableVector",init<int>())
        .def("__len__",&TV::size)
        .def("__getitem__",&get_item<TV,N,T>)
        .def("__setitem__",&set_item<TV,N,T>)
        .def("__setitem__",&set_item<TV,N,T>)
        .def(self_ns::str(self))
        ;


/*
      def("asin", (T(*)(const T&))&asin);
      def("acos", (T(*)(const T&))&acos);
      def("atan", (T(*)(const T&))&atan);
    */

}

void calculus_submodule() 
{
    export_taylor_variable();
}


