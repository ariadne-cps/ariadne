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

TaylorVariable mul_cosy(const TaylorVariable&, const TaylorVariable&);
TaylorVariable mul_rounded(const TaylorVariable&, const TaylorVariable&);
TaylorVariable mul_ivl(const TaylorVariable&, const TaylorVariable&);
    
    
void read(MultiIndex& j, const boost::python::object& obj) {
    array<uint> ary;
    read_tuple_array(ary,obj);
    j=MultiIndex(ary.size(),ary.begin());
}


void read(SparseDifferential<Float>& sd, const boost::python::object& obj) {
    Float c;
    MultiIndex j(0);
    boost::python::dict dct=extract<boost::python::dict>(obj);
    boost::python::list lst=dct.items();
    for(uint i=0; i!=len(lst); ++i) {
        boost::python::tuple tup=extract<boost::python::tuple>(lst[i]);
        read(j,tup[0]);
        read(c,tup[1]);
        if(i==0) { 
            sd=SparseDifferential<Float>(j.size());
        }
        sd[j]=c;
    }
}



void read(TaylorVariable& tv, const boost::python::object& obj1, const boost::python::object& obj2) {
    SparseDifferential<Float> sd;
    Interval e;
    read(sd,obj1);
    read(e,obj2); 
    if(boost::python::extract<Float>(obj2).check() && e.u>=0) { e.l=-e.u; } 
    else if(e.l!=-e.u) { std::cerr<<"WARNING: TaylorVariable error is not a symmetric interval.\n"; }
    tv=TaylorVariable(sd,e);
}
   

void read(TaylorVariable& tv, const boost::python::object& obj1, const boost::python::object& obj2,
          const boost::python::object& obj3, const boost::python::object& obj4) 
{
    uint as;
    uint deg;
    array<Float> dat;
    Interval err;
    read(as,obj1);
    read(deg,obj2);
    read_list_array(dat,obj3);
    read(err,obj4);
    const Float* ptr=dat.begin();
    if(boost::python::extract<Float>(obj4).check() && err.u>=0) { err.l=-err.u; } 
    else if(err.l!=-err.u) { std::cerr<<"WARNING: TaylorVariable error is not a symmetric interval.\n"; }
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
        SparseDifferential<Float> sd;
        Interval e=0;
        read(sd,obj);
        tv=TaylorVariable(sd,e);
    } 
}
   
TaylorVariable*
make_taylor_variable(const uint& as, const uint& deg, const boost::python::object& aobj, const boost::python::object& eobj) 
{
    array<Float> dat;
    Interval err;
    read_list_array(dat,aobj);
    const Float* ptr=dat.begin();
    read(err,eobj);
    if(boost::python::extract<Float>(eobj).check() && err.u>=0) { err.l=-err.u; } 
    else if(err.l!=-err.u) { std::cerr<<"WARNING: TaylorVariable error is not a symmetric interval.\n"; }
    return new TaylorVariable(as,deg,ptr,err);
}

boost::python::list
make_taylor_variables(const Vector<Interval>& x) 
{
    boost::python::list result;
    for(uint i=0; i!=x.size(); ++i) {
        result.append(TaylorVariable::variable(x.size(),midpoint(x[i]),i));
    }
    return result;
}

   


boost::python::tuple
python_split(const TaylorVariable& x, uint j)
{
    std::pair<TaylorVariable,TaylorVariable> res=split(x,j);
    return boost::python::make_tuple(res.first,res.second);
}

template<class C, class I, class X> inline 
void set_item(C& c, const I& i, const X& x) {
    c[i]=x;
}


template<class C, class I> inline 
typename C::ValueType 
get_item(const C& c, const I& i) {
    return c[i];
}

template<> std::string __repr__(const TaylorVariable& tv) {
    std::stringstream ss;
    ss << "TV(" << tv.expansion() << "," << tv.error() << ")";
    return ss.str();
} 

} // namespace Ariadne

void export_taylor_variable() 
{
    typedef double D;
    typedef Float R;
    typedef Interval I;
    typedef MultiIndex A;
    typedef Vector<Float> V;
    typedef TaylorVariable T;


    class_<T> taylor_variable_class("TaylorVariable");
    taylor_variable_class.def("__init__", make_constructor(&make<TaylorVariable>) );
    taylor_variable_class.def("__init__", make_constructor(&make2<TaylorVariable>) );
    taylor_variable_class.def("__init__", make_constructor(&make_taylor_variable) );
    taylor_variable_class.def( init< uint >());
    taylor_variable_class.def("error", (const I&(T::*)()const) &T::error, return_value_policy<copy_const_reference>());
    taylor_variable_class.def("clean", &TaylorVariable::clean);
    taylor_variable_class.def("domain", &TaylorVariable::domain);
    taylor_variable_class.def("range", &TaylorVariable::range);
    taylor_variable_class.def("__getitem__", &get_item<T,A>);
    taylor_variable_class.def("__setitem__",&set_item<T,A,D>);
    taylor_variable_class.def("__setitem__",&set_item<T,A,R>);
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
    taylor_variable_class.def(self_ns::str(self));
    taylor_variable_class.def("truncate", &TaylorVariable::truncate,return_value_policy<reference_existing_object>());
    taylor_variable_class.def("sweep", &TaylorVariable::sweep,return_value_policy<reference_existing_object>());
    taylor_variable_class.def("evaluate", (Interval(TaylorVariable::*)(const Vector<Float>&)const) &TaylorVariable::evaluate);
    taylor_variable_class.def("evaluate", (Interval(TaylorVariable::*)(const Vector<Interval>&)const) &TaylorVariable::evaluate);

    def("scale", (T(*)(const T&,const Interval&)) &scale);
    def("split", &python_split);


    taylor_variable_class.def("__repr__",&__repr__<T>);


    taylor_variable_class.def("constant",(T(*)(uint, const R&))&T::constant);
    taylor_variable_class.def("variable",(T(*)(uint, const R&, uint))&T::variable);
    taylor_variable_class.def("variables",&make_taylor_variables);
  
    taylor_variable_class.staticmethod("constant");
    taylor_variable_class.staticmethod("variable");
    taylor_variable_class.staticmethod("variables");

    def("mul_cosy", (T(*)(const T&, const T&)) &mul_cosy);
    def("mul_rounded", (T(*)(const T&, const T&)) &mul_rounded);
    def("mul_ivl", (T(*)(const T&, const T&)) &mul_ivl);

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


