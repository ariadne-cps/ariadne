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

#include "expression_interface.h"
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

TaylorExpression midpoint(const TaylorExpression& x) {
    TaylorExpression r=x;
    r.set_error(0.0);
    return r;
}

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


void read(TaylorExpression& te, const boost::python::object& obj);

void read(TaylorExpression& te, const boost::python::tuple& tup) {
    if(len(tup)==1) {
        read(te,tup[0]);
    }else if(len(tup)==3) {
        read(te,tup[0],tup[1],tup[2]);
    }
}


void read(TaylorExpression& te, const boost::python::object& obj) {
    if(check(extract<boost::python::tuple>(obj))) {
         read(te,extract<boost::python::tuple>(obj));
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



template<> std::string __repr__(const TaylorExpression& te) {
    std::stringstream ss;
    ss << "TaylorExpression(" << te.expansion() << "," << te.error() << ")";
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
    Float r;
    if(pi.expansion().size()==0) { r=0.0; } else { r=radius(pi.begin()->data()); }
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

std::string __str__(const TaylorExpression& te) {
    std::stringstream ss;
    Polynomial<Interval> p=te.polynomial();
    ss<<"TaylorExpression"<<__tbox_str__(te.domain())<<"( "<<__tpoly_str__(p)<<" )";
    ss<<"( m=" << te.model()<<" )";
    return ss.str();
}


TaylorFunction __getslice__(const TaylorFunction& tf, int start, int stop) {
    return TaylorFunction(tf.domain(),Vector<TaylorModel>(project(tf.models(),range(start,stop))));
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
    typedef TaylorExpression TE;
    typedef TaylorFunction TF;
    typedef Polynomial<Float> RP;
    typedef Polynomial<Interval> IP;

    class_<TE> taylor_expression_class("TaylorExpression");
    taylor_expression_class.def("__init__", make_constructor(&make<TE>) );
    taylor_expression_class.def("__init__", make_constructor(&make3<TE>) );
    taylor_expression_class.def( init< IV >());
    taylor_expression_class.def( init< IV, RP >());
    taylor_expression_class.def("error", (const R&(TE::*)()const) &TE::error, return_value_policy<copy_const_reference>());
    taylor_expression_class.def("argument_size", &TE::argument_size);
    taylor_expression_class.def("domain", &TE::domain, return_value_policy<copy_const_reference>());
    taylor_expression_class.def("centre", &TE::centre);
    taylor_expression_class.def("range", &TE::range);
    taylor_expression_class.def("__getitem__", &get_item<TE,A,R>);
    taylor_expression_class.def("__setitem__",&set_item<TE,A,D>);
    taylor_expression_class.def(-self);
    taylor_expression_class.def(self+self);
    taylor_expression_class.def(self-self);
    taylor_expression_class.def(self*self);
    taylor_expression_class.def(self/self);
    taylor_expression_class.def(self+R());
    taylor_expression_class.def(self-R());
    taylor_expression_class.def(self*R());
    taylor_expression_class.def(self/R());
    taylor_expression_class.def(R()+self);
    taylor_expression_class.def(R()-self);
    taylor_expression_class.def(R()*self);
    taylor_expression_class.def(R()/self);
    taylor_expression_class.def(self+=R());
    taylor_expression_class.def(self-=R());
    taylor_expression_class.def(self*=R());
    taylor_expression_class.def(self/=R());
    taylor_expression_class.def(self+I());
    taylor_expression_class.def(self-I());
    taylor_expression_class.def(self*I());
    taylor_expression_class.def(self/I());
    taylor_expression_class.def(self+=I());
    taylor_expression_class.def(self-=I());
    taylor_expression_class.def(self*=I());
    taylor_expression_class.def(self/=I());
    taylor_expression_class.def(self+=self);
    taylor_expression_class.def(self-=self);
    taylor_expression_class.def(self>R());
    taylor_expression_class.def(self>self);
    taylor_expression_class.def(self<self);
    //taylor_expression_class.def(self_ns::str(self));
    taylor_expression_class.def("__str__",(std::string(*)(const TE&)) &__str__);
    taylor_expression_class.def("value", (const Float&(TE::*)()const) &TE::value,return_value_policy<copy_const_reference>());
    taylor_expression_class.def("truncate", (TE&(TE::*)(uint)) &TE::truncate,return_value_policy<reference_existing_object>());
    taylor_expression_class.def("sweep", (TE&(TE::*)(double))&TE::sweep,return_value_policy<reference_existing_object>());
    taylor_expression_class.def("truncate", (TE&(TE::*)()) &TE::truncate,return_value_policy<reference_existing_object>());
    taylor_expression_class.def("sweep", (TE&(TE::*)())&TE::sweep,return_value_policy<reference_existing_object>());
    taylor_expression_class.def("clobber", (TE&(TE::*)()) &TE::clobber,return_value_policy<reference_existing_object>());
    taylor_expression_class.def("clean", (TE&(TE::*)()) &TE::clean,return_value_policy<reference_existing_object>());
    taylor_expression_class.def("set_maximum_degree",&TE::set_maximum_degree);
    taylor_expression_class.def("set_sweep_threshold",&TE::set_sweep_threshold);
    taylor_expression_class.def("maximum_degree",&TE::maximum_degree);
    taylor_expression_class.def("sweep_threshold",&TE::sweep_threshold);
    //taylor_expression_class.def("set_default_maximum_degree",&TE::set_default_maximum_degree);
    //taylor_expression_class.def("set_default_sweep_threshold",&TE::set_default_sweep_threshold);
    //taylor_expression_class.def("default_maximum_degree",&TE::default_maximum_degree);
    //taylor_expression_class.def("default_sweep_threshold",&TE::default_sweep_threshold);
    taylor_expression_class.def("__call__", (Interval(TE::*)(const Vector<Float>&)const) &TE::evaluate);
    taylor_expression_class.def("__call__", (Interval(TE::*)(const Vector<Interval>&)const) &TE::evaluate);
    taylor_expression_class.def("evaluate", (Interval(TE::*)(const Vector<Float>&)const) &TE::evaluate);
    taylor_expression_class.def("evaluate", (Interval(TE::*)(const Vector<Interval>&)const) &TE::evaluate);
    taylor_expression_class.def("polynomial", (Polynomial<Interval>(TE::*)()const) &TE::polynomial);
    //taylor_expression_class.staticmethod("set_default_maximum_degree");
    //taylor_expression_class.staticmethod("set_default_sweep_threshold");
    //taylor_expression_class.staticmethod("default_maximum_degree");
    //taylor_expression_class.staticmethod("default_sweep_threshold");
    def("split", &python_split);


    taylor_expression_class.def("__repr__",&__repr__<TE>);

    taylor_expression_class.def("constant",(TE(*)(const IV&, const R&))&TE::constant);
    taylor_expression_class.def("variable",(TE(*)(const IV&, uint))&TE::variable);
    taylor_expression_class.def("variables",(TF(*)(const IV&))&TE::variables);

    taylor_expression_class.staticmethod("constant");
    taylor_expression_class.staticmethod("variable");
    taylor_expression_class.staticmethod("variables");

    def("midpoint",(TE(*)(const TE&)) &midpoint);
    def("refines",(bool(*)(const TE&,const TE&)) &refines);
    def("evaluate",(I(*)(const TE&,const IV&)) &evaluate);
    def("derivative",(TE(*)(const TE&,N)) &derivative);
    def("antiderivative",(TE(*)(const TE&,N)) &antiderivative);

    def("embed",(TE(*)(const TE&,const Interval&)) &embed);
    def("embed",(TE(*)(const TE&,const IV&)) &embed);
    def("embed",(TE(*)(const IV&,const TE&)) &embed);
    //def("/*compose*/",(TE(*)(const TE&,const TE&)) &compose);
    //def("compose",(T(*)(const T&,const TE&)) &compose);

    def("max",(TE(*)(const TE&,const TE&))&max);
    def("min",(TE(*)(const TE&,const TE&))&min);
    def("abs",(TE(*)(const TE&))&abs);

    def("neg",(TE(*)(const TE&))&neg);
    def("rec",(TE(*)(const TE&))&rec);
    def("sqr",(TE(*)(const TE&))&sqr);
    def("pow",(TE(*)(const TE&, int))&pow);

    def("sqrt", (TE(*)(const TE&))&sqrt);
    def("exp", (TE(*)(const TE&))&exp);
    def("log", (TE(*)(const TE&))&log);
    def("sin", (TE(*)(const TE&))&sin);
    def("cos", (TE(*)(const TE&))&cos);
    def("tan", (TE(*)(const TE&))&tan);

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
    typedef TaylorExpression TE;
    typedef TaylorFunction TF;
    typedef ExpressionInterface E;
    typedef FunctionInterface F;

    class_<TF> taylor_function_class("TaylorFunction");
    taylor_function_class.def( init< IV >());
    taylor_function_class.def( init< IV,const F& >());
    taylor_function_class.def("result_size", &TaylorFunction::result_size);
    taylor_function_class.def("argument_size", &TaylorFunction::argument_size);
    taylor_function_class.def("domain", &TaylorFunction::domain, return_value_policy<copy_const_reference>());
    taylor_function_class.def("centre", &TaylorFunction::centre);
    taylor_function_class.def("range", &TaylorFunction::range);
    taylor_function_class.def("__getslice__", &__getslice__);
    taylor_function_class.def("__getitem__", &get_item<TF,N,TE>);
    taylor_function_class.def("__setitem__",&set_item<TF,N,TE>);
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
    taylor_function_class.def("clobber", (TF&(TF::*)()) &TF::clobber,return_value_policy<reference_existing_object>());
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
    taylor_function_class.def("polynomial", (Vector< Polynomial<Interval> >(TF::*)()const) &TF::polynomial);
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
    def("join", (TF(*)(const TF&,const TE&)) &join);
    def("join", (TF(*)(const TE&,const TE&)) &join);

    def("combine", (TF(*)(const TF&,const TE&)) &combine);
    def("combine", (TF(*)(const TF&,const TF&)) &combine);

    def("embed",(TF(*)(const TF&,const Interval&)) &embed);
    def("embed",(TF(*)(const TF&,const IV&)) &embed);
    def("embed",(TF(*)(const IV&,const TF&)) &embed);

    def("evaluate",(IV(TF::*)(const RV&)const) &TF::evaluate);
    def("evaluate",(IV(TF::*)(const IV&)const) &TF::evaluate);
    def("compose",(TE(*)(const TE&,const TF&)) &compose);
    def("compose",(TF(*)(const TF&,const TF&)) &compose);
    def("compose",(TE(*)(const E&,const TF&)) &compose);
    def("compose",(TF(*)(const F&,const TF&)) &compose);
    def("implicit",(TE(*)(const TE&)) &implicit);
    def("implicit",(TF(*)(const TF&)) &implicit);
    def("antiderivative",(TF(*)(const TF&,N)) &antiderivative);
    def("flow",(TF(*)(const TF&,const IV&,const I&, N)) &flow);

    def("unchecked_compose",(TE(*)(const TE&,const TF&)) &unchecked_compose);
    def("unchecked_compose",(TF(*)(const TF&,const TF&)) &unchecked_compose);
    def("unchecked_flow",(TF(*)(const TF&,const IV&,const I&, N)) &unchecked_flow);



}

void calculus_submodule()
{
    export_taylor_variable();
    export_taylor_function();
}


