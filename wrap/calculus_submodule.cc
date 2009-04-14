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

#include "function_interface.h"
#include "polynomial.h"
#include "taylor_expression.h"
#include "taylor_function.h"

#include <boost/python.hpp>
using namespace boost::python;

using namespace Ariadne;

#include "utilities.h"

namespace Ariadne {

TaylorExpression max(const TaylorExpression& x, const TaylorExpression& y);
TaylorExpression min(const TaylorExpression& x, const TaylorExpression& y);

void read(MultiIndex& j, const boost::python::object& obj) {
    array<int> ary;
    read_tuple_array(ary,obj);
    j=MultiIndex(ary.size(),ary.begin());
}


template<class K, class V>
void read(std::map<K,V>& m, const boost::python::object& obj) {
    V c;
    K j;
    boost::python::dict dct=extract<boost::python::dict>(obj);
    boost::python::list lst=dct.items();
    for(int i=0; i!=len(lst); ++i) {
        boost::python::tuple tup=extract<boost::python::tuple>(lst[i]);
        read(j,tup[0]);
        read(c,tup[1]);
        if(i==0) {
            m=std::map<K,V>();
        }
        m[j]=c;
    }
}


template<class X>
void read(Vector<X>& v, const boost::python::object& obj);


void read(TaylorModel& t, const boost::python::object& obj1, const boost::python::object& obj2) {
    std::map<MultiIndex,Float> m;
    Float e;
    read(m,obj1);
    read(e,obj2);
    ARIADNE_ASSERT(e>=0);
    t=TaylorModel(m,e);
}


void read(TaylorExpression& t, const boost::python::object& obj1, const boost::python::object& obj2, const boost::python::object& obj3) {
    Vector<Interval> d;
    std::map<MultiIndex,Float> m;
    Float e;
    read(d,obj1);
    read(m,obj2);
    read(e,obj3);
    ARIADNE_ASSERT(e>=0);
    t=TaylorExpression(d,m,e);
}


void read(TaylorExpression& tv, const boost::python::object& obj);

void read(TaylorExpression& tv, const boost::python::tuple& tup) {
    if(len(tup)==1) {
        read(tv,tup[0]);
    }else if(len(tup)==3) {
        read(tv,tup[0],tup[1],tup[2]);
    }
}


void read(TaylorExpression& tv, const boost::python::object& obj) {
    if(check(extract<boost::python::tuple>(obj))) {
         read(tv,extract<boost::python::tuple>(obj));
    }
}

Vector<TaylorExpression>
make_taylor_variables(const Vector<Interval>& x)
{
    return TaylorExpression::variables(x);
}




boost::python::tuple
python_split(const TaylorExpression& x, uint j)
{
    std::pair<TaylorExpression,TaylorExpression> res=split(x,j);
    return boost::python::make_tuple(res.first,res.second);
}



template<> std::string __repr__(const TaylorExpression& tv) {
    std::stringstream ss;
    ss << "TaylorExpression(" << tv.expansion() << "," << tv.error() << ")";
    return ss.str();
}

std::string __tbox_str__(const Vector<Interval>& d) {
    std::stringstream ss;
    for(uint j=0; j!=d.size(); ++j) {
        ss<<(j==0?"[":"x")<<d[j];
    }
    ss<<"]";
    return ss.str();
}


std::string __tpoly_str__(const Polynomial<Interval>& pi) {
    std::stringstream ss;
    bool first=true;
    Float r=radius(pi.begin()->data());
    for(Polynomial<Interval>::const_iterator iter=pi.begin(); iter!=pi.end(); ++iter) {
        MultiIndex a=iter->key();
        Float v=midpoint(iter->data());
        if(abs(v)<1e-15) { r+=abs(v); v=0; }
        if(v!=0) {
            if(v>0 && !first) { ss<<"+"; }
            first=false;
            if(v<0) { ss<<"-"; }
            if(abs(v)!=1) { ss<<abs(v); }
            for(uint j=0; j!=a.size(); ++j) {
                if(a[j]!=0) { ss<<"x"<<j; if(a[j]!=1) { ss<<"^"<<int(a[j]); } }
            }
        }
    }
    if(r>0) { ss<<"+/-"<<r; }
    return ss.str();
}

std::string __str__(const TaylorExpression& tv) {
    std::stringstream ss;
    Polynomial<Interval> p=TaylorFunction(tv.domain(),Vector<TaylorModel>(1,tv.model())).polynomial()[0];
    ss<<"TaylorExpression"<<__tbox_str__(tv.domain())<<"( "<<__tpoly_str__(p)<<" )";
    return ss.str();
}



std::string __str__(const TaylorFunction& tf) {
    std::stringstream ss;
    Vector< Polynomial<Interval> > p=tf.polynomial();
    ss<<"TaylorFunction"<<__tbox_str__(tf.domain());
    for(uint i=0; i!=p.size(); ++i) {
        ss<<(i==0?"( ":"; ")<<__tpoly_str__(p[i]);
    }
    ss<<" )";
    return ss.str();

    return ss.str();
}




TaylorExpression neg(const TaylorExpression&);
TaylorExpression rec(const TaylorExpression&);
TaylorExpression sqr(const TaylorExpression&);
TaylorExpression pow(const TaylorExpression&, int);
TaylorExpression sqrt(const TaylorExpression&);
TaylorExpression exp(const TaylorExpression&);
TaylorExpression log(const TaylorExpression&);
TaylorExpression sin(const TaylorExpression&);
TaylorExpression cos(const TaylorExpression&);
TaylorExpression tan(const TaylorExpression&);

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
    typedef TaylorExpression TV;
    typedef TaylorFunction TF;

    class_<TV> taylor_variable_class("TaylorExpression");
    taylor_variable_class.def("__init__", make_constructor(&make<TV>) );
    taylor_variable_class.def("__init__", make_constructor(&make3<TV>) );
    taylor_variable_class.def( init< IV >());
    taylor_variable_class.def( init< IV >());
    taylor_variable_class.def( init< TV >());
    taylor_variable_class.def("error", (const R&(TV::*)()const) &TV::error, return_value_policy<copy_const_reference>());
    taylor_variable_class.def("domain", &TV::domain, return_value_policy<copy_const_reference>());
    taylor_variable_class.def("range", &TV::range);
    taylor_variable_class.def("__getitem__", &get_item<TV,A,R>);
    taylor_variable_class.def("__setitem__",&set_item<TV,A,D>);
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
    taylor_variable_class.def(self>R());
    taylor_variable_class.def(self>self);
    taylor_variable_class.def(self<self);
    //taylor_variable_class.def(self_ns::str(self));
    taylor_variable_class.def("__str__",(std::string(*)(const TV&)) &__str__);
    taylor_variable_class.def("truncate", (TV&(TV::*)(uint)) &TV::truncate,return_value_policy<reference_existing_object>());
    taylor_variable_class.def("sweep", (TV&(TV::*)(double))&TV::sweep,return_value_policy<reference_existing_object>());
    taylor_variable_class.def("truncate", (TV&(TV::*)()) &TV::truncate,return_value_policy<reference_existing_object>());
    taylor_variable_class.def("sweep", (TV&(TV::*)())&TV::sweep,return_value_policy<reference_existing_object>());
    taylor_variable_class.def("clean", (TV&(TV::*)()) &TV::clean,return_value_policy<reference_existing_object>());
    taylor_variable_class.def("set_maximum_degree",&TV::set_maximum_degree);
    taylor_variable_class.def("set_sweep_threshold",&TV::set_sweep_threshold);
    taylor_variable_class.def("maximum_degree",&TV::maximum_degree);
    taylor_variable_class.def("sweep_threshold",&TV::sweep_threshold);
    //taylor_variable_class.def("set_default_maximum_degree",&TV::set_default_maximum_degree);
    //taylor_variable_class.def("set_default_sweep_threshold",&TV::set_default_sweep_threshold);
    //taylor_variable_class.def("default_maximum_degree",&TV::default_maximum_degree);
    //taylor_variable_class.def("default_sweep_threshold",&TV::default_sweep_threshold);
    taylor_variable_class.def("__call__", (Interval(TV::*)(const Vector<Float>&)const) &TV::evaluate);
    taylor_variable_class.def("__call__", (Interval(TV::*)(const Vector<Interval>&)const) &TV::evaluate);
    taylor_variable_class.def("evaluate", (Interval(TV::*)(const Vector<Float>&)const) &TV::evaluate);
    taylor_variable_class.def("evaluate", (Interval(TV::*)(const Vector<Interval>&)const) &TV::evaluate);
    //taylor_variable_class.staticmethod("set_default_maximum_degree");
    //taylor_variable_class.staticmethod("set_default_sweep_threshold");
    //taylor_variable_class.staticmethod("default_maximum_degree");
    //taylor_variable_class.staticmethod("default_sweep_threshold");
    def("split", &python_split);


    taylor_variable_class.def("__repr__",&__repr__<TV>);

    taylor_variable_class.def("constant",(TV(*)(const IV&, const R&))&TV::constant);
    taylor_variable_class.def("variable",(TV(*)(const IV&, uint))&TV::variable);
    taylor_variable_class.def("variables",(TF(*)(const IV&))&TV::variables);

    taylor_variable_class.staticmethod("constant");
    taylor_variable_class.staticmethod("variable");
    taylor_variable_class.staticmethod("variables");

    def("evaluate",(I(*)(const TV&,const IV&)) &evaluate);

    //def("compose",(TV(*)(const TV&,const TV&)) &compose);
    //def("compose",(T(*)(const T&,const TV&)) &compose);

    def("max",(TV(*)(const TV&,const TV&))&max);
    def("min",(TV(*)(const TV&,const TV&))&min);
    def("abs",(TV(*)(const TV&))&abs);

    def("neg",(TV(*)(const TV&))&neg);
    def("rec",(TV(*)(const TV&))&rec);
    def("sqr",(TV(*)(const TV&))&sqr);
    def("pow",(TV(*)(const TV&, int))&pow);

    def("sqrt", (TV(*)(const TV&))&sqrt);
    def("exp", (TV(*)(const TV&))&exp);
    def("log", (TV(*)(const TV&))&log);
    def("sin", (TV(*)(const TV&))&sin);
    def("cos", (TV(*)(const TV&))&cos);
    def("tan", (TV(*)(const TV&))&tan);

}

void export_taylor_function()
{
    typedef uint N;
    typedef double D;
    typedef Float R;
    typedef Interval I;
    typedef MultiIndex A;
    typedef Vector<Float> RV;
    typedef Vector<Interval> IV;
    typedef Matrix<Float> RMx;
    typedef TaylorExpression TV;
    typedef TaylorFunction TF;
    typedef FunctionInterface F;

    class_<TF> taylor_function_class("TaylorFunction");
    taylor_function_class.def( init< IV >());
    taylor_function_class.def( init< IV,const F& >());
    taylor_function_class.def("result_size", &TaylorFunction::result_size);
    taylor_function_class.def("argument_size", &TaylorFunction::argument_size);
    taylor_function_class.def("domain", &TaylorFunction::domain, return_value_policy<copy_const_reference>());
    taylor_function_class.def("range", &TaylorFunction::range);
    taylor_function_class.def("__getitem__", &get_item<TF,N,TV>);
    taylor_function_class.def("__setitem__",&set_item<TF,N,TV>);
    taylor_function_class.def(-self);
    taylor_function_class.def(self+self);
    taylor_function_class.def(self-self);
    taylor_function_class.def(self+RV());
    taylor_function_class.def(self-RV());
    taylor_function_class.def(self*R());
    taylor_function_class.def(self/R());
    //taylor_function_class.def(RV()+self);
    //taylor_function_class.def(RV()-self);
    //taylor_function_class.def(R()*self);
    //taylor_function_class.def(R()/self);
    taylor_function_class.def(self+=RV());
    taylor_function_class.def(self-=RV());
    taylor_function_class.def(self*=R());
    taylor_function_class.def(self/=R());
    taylor_function_class.def(self+=self);
    taylor_function_class.def(self-=self);
    taylor_function_class.def("__str__",(std::string(*)(const TF&)) &__str__);
    //taylor_function_class.def(self_ns::str(self));
    //taylor_function_class.def("truncate", (TaylorFunction&(TaylorFunction::*)(uint)) &TaylorFunction::truncate,return_value_policy<reference_existing_object>());
    //taylor_function_class.def("sweep", (TaylorFunction&(TaylorFunction::*)(double))&TaylorFunction::sweep,return_value_policy<reference_existing_object>());
    //taylor_function_class.def("truncate", (TaylorFunction&(TaylorFunction::*)()) &TaylorFunction::truncate,return_value_policy<reference_existing_object>());
    //taylor_function_class.def("sweep", (TaylorFunction&(TaylorFunction::*)())&TaylorFunction::sweep,return_value_policy<reference_existing_object>());
    //taylor_function_class.def("clean", (TaylorFunction&(TaylorFunction::*)()) &TaylorFunction::clean,return_value_policy<reference_existing_object>());
    //taylor_function_class.def("set_maximum_degree",&TaylorFunction::set_maximum_degree);
    //taylor_function_class.def("set_sweep_threshold",&TaylorFunction::set_sweep_threshold);
    //taylor_function_class.def("maximum_degree",&TaylorFunction::maximum_degree);
    //taylor_function_class.def("sweep_threshold",&TaylorFunction::sweep_threshold);
    //taylor_function_class.def("set_default_maximum_degree",&TaylorFunction::set_default_maximum_degree);
    //taylor_function_class.def("set_default_sweep_threshold",&TaylorFunction::set_default_sweep_threshold);
    //taylor_function_class.def("default_maximum_degree",&TaylorFunction::default_maximum_degree);
    //taylor_function_class.def("default_sweep_threshold",&TaylorFunction::default_sweep_threshold);
    taylor_function_class.def("__call__", (Vector<Interval>(TaylorFunction::*)(const Vector<Float>&)const) &TaylorFunction::evaluate);
    taylor_function_class.def("__call__", (Vector<Interval>(TaylorFunction::*)(const Vector<Interval>&)const) &TaylorFunction::evaluate);
    taylor_function_class.def("evaluate", (Vector<Interval>(TaylorFunction::*)(const Vector<Float>&)const) &TaylorFunction::evaluate);
    taylor_function_class.def("evaluate", (Vector<Interval>(TaylorFunction::*)(const Vector<Interval>&)const) &TaylorFunction::evaluate);
    //taylor_function_class.staticmethod("set_default_maximum_degree");
    //taylor_function_class.staticmethod("set_default_sweep_threshold");
    //taylor_function_class.staticmethod("default_maximum_degree");
    //taylor_function_class.staticmethod("default_sweep_threshold");

    taylor_function_class.def("__repr__",&__repr__<TF>);

    taylor_function_class.def("constant",(TF(*)(const IV&, const RV&))&TF::constant);
    taylor_function_class.def("constant",(TF(*)(const IV&, const IV&))&TF::constant);
    taylor_function_class.def("identity",(TF(*)(const IV&))&TF::identity);

    taylor_function_class.staticmethod("constant");
    taylor_function_class.staticmethod("identity");
 
    def("join", (TF(*)(const TF&,const TF&)) &join);
    def("join", (TF(*)(const TF&,const TV&)) &join);

    def("combine", (TF(*)(const TF&,const TV&)) &combine);
    def("combine", (TF(*)(const TF&,const TF&)) &combine);

    def("evaluate",(IV(TF::*)(const RV&)const) &TF::evaluate);
    def("evaluate",(IV(TF::*)(const IV&)const) &TF::evaluate);
    def("compose",(TV(*)(const TV&,const TF&)) &compose);
    def("compose",(TF(*)(const TF&,const TF&)) &compose);
    def("compose",(TF(*)(const F&,const TF&)) &compose);
    def("implicit",(TV(*)(const TV&)) &implicit);
    def("implicit",(TF(*)(const TF&)) &implicit);
    def("antiderivative",(TF(*)(const TF&,N)) &antiderivative);
    def("flow",(TF(*)(const TF&,const IV&,const I&, N)) &flow);
    def("unchecked_flow",(TF(*)(const TF&,const IV&,const I&, N)) &unchecked_flow);



}

void calculus_submodule()
{
    export_taylor_variable();
    export_taylor_function();
}


