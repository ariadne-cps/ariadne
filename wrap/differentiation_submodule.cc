/***************************************************************************
 *            differentiation_submodule.cc
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
 
#include "array.h"
#include "numeric.h"
#include "taylor_model.h"
#include "differential.h"

#include <boost/python.hpp>
using namespace boost::python;

using namespace Ariadne;

template<class X> void read_array(array<X>&, const boost::python::object& obj) { }
inline uint compute_polynomial_data_size(uint rs, uint as, uint d) { return rs*Ariadne::bin(d+as,as); }

template<class DIFF>
DIFF*
make_differential(const uint& as, const uint& d, const boost::python::object& obj) 
{
    typedef typename DIFF::scalar_type X;
    DIFF* result=new DIFF(as,d);
    array<X> data;
    read_array(data,obj);
    assert(data.size()==compute_polynomial_data_size(1u,as,d));
    MultiIndex i(as);
    const X* ptr=data.begin();
    while(i.degree()<=d) {
        result[i]=*ptr; ++i; ++ptr;
    }
    return result;
}


template<class DIFF>
boost::python::list
make_differential_variables(const uint& d, const Vector<Interval>& x) 
{
    boost::python::list result;
    for(uint i=0; i!=x.size(); ++i) {
        result.append(DIFF::variable(x.size(),d,midpoint(x[i]),i));
    }
    return result;
}


template<class DIFF>
Vector<DIFF>*
make_differential_vector(const uint& rs, const uint& as, const uint& d, const boost::python::object& obj) 
{
    typedef typename DIFF::scalar_type X;
    array<X> data;
    read_array(data,obj);
    ARIADNE_ASSERT(data.size()==compute_polynomial_data_size(rs,as,d));
    Vector<DIFF>* result=new Vector<DIFF>(rs,DIFF(as,d,data.begin()));
    return result;
}


template<class C, class I, class X> inline 
X get_item(const C& c, const I& i) { return c[i]; }

template<class C, class I, class J, class X> inline 
X matrix_get_item(const C& c, const I& i, const J& j) { return c[i][j]; }

template<class C, class I, class X> inline 
void set_item(C& c, const I& i, const X& x) { c[i]=x; }

template<class C, class I, class J, class X> inline 
void matrix_set_item(C& c, const I& i, const J& j, const X& x) { c[i][j]=x; }




template<class DIFF>
void export_differential(const char* name) 
{
    typedef typename DIFF::scalar_type X;
    typedef Vector<X> V;
    typedef Series<X> S;
    typedef DIFF D;
    typedef Vector<D> DV;


    class_<D> differential_class(name);
    //differential_class.def("__init__", make_constructor(&make_differential<X>) );
    differential_class.def( init< uint, uint >());
    differential_class.def("value", (const X&(D::*)()const) &D::value, return_value_policy<copy_const_reference>());
    differential_class.def("__getitem__", &get_item<D,MultiIndex,X>);
    differential_class.def("__setitem__",&set_item<D,MultiIndex,double>);
    differential_class.def("__setitem__",&set_item<D,MultiIndex,X>);
    differential_class.def(-self);
    differential_class.def(self+self);
    differential_class.def(self-self);
    differential_class.def(self*self);
    differential_class.def(self/self);
    differential_class.def(self+X());
    differential_class.def(self-X());
    differential_class.def(self+=X());
    differential_class.def(self-=X());
    differential_class.def(self*=X());
    differential_class.def(self/=X());
    differential_class.def(self_ns::str(self));
  
    differential_class.def("constant",(D(*)(uint, ushort, const X&))&D::constant);
    differential_class.def("variable",(D(*)(uint, ushort, const X&, uint))&D::variable);
    differential_class.def("variables",&make_differential_variables<D>);

    differential_class.staticmethod("constant");
    differential_class.staticmethod("variable");
    differential_class.staticmethod("variables");

    def("max",(D(*)(const D&,const D&))&max<X>);
    def("min",(D(*)(const D&,const D&))&min<X>);
    def("abs",(D(*)(const D&))&abs<X>);
    def("pos",(D(*)(const D&))&pos<X>);
    def("neg",(D(*)(const D&))&neg<X>);
    def("rec",(D(*)(const D&))&rec<X>);
    def("pow",(D(*)(const D&, int))&pow<X>);

    def("sqrt", (D(*)(const D&))&sqrt<X>);
    def("exp", (D(*)(const D&))&exp<X>);
    def("log", (D(*)(const D&))&log<X>);
/*
    def("sin", (D(*)(const D&))&sin<X>);
    def("cos", (D(*)(const D&))&cos<X>);
    def("tan", (D(*)(const D&))&tan<X>);
    def("asin", (D(*)(const D&))&asin<X>);
    def("acos", (D(*)(const D&))&acos<X>);
    def("atan", (D(*)(const D&))&atan<X>);
*/
}

template<class DIFF> 
void
export_differential_vector(const char* name)
{
    typedef typename DIFF::scalar_type X;
    typedef Vector<X> V;
    typedef Series<X> S;
    typedef DIFF D;
    typedef Vector<D> DV;

    class_<DV> differential_vector_class(name);
    differential_vector_class.def("__init__", make_constructor(&make_differential_vector<D>) );
    differential_vector_class.def( init< uint, uint, uint >());
    differential_vector_class.def("__getitem__", &matrix_get_item<DV,int,MultiIndex,X>);
    differential_vector_class.def("__setitem__",&matrix_set_item<DV,int,MultiIndex,double>);
    differential_vector_class.def("__getitem__", &get_item<DV,int,D>);
    differential_vector_class.def("__setitem__",&set_item<DV,int,X>);
    differential_vector_class.def("__setitem__",&set_item<DV,int,D>);
    differential_vector_class.def(-self);
    differential_vector_class.def(self+self);
    differential_vector_class.def(self-self);
    differential_vector_class.def(self+V());
    differential_vector_class.def(self-V());
    //differential_vector_class.def(self*X());
    differential_vector_class.def(self+=V());
    differential_vector_class.def(self-=V());
    differential_vector_class.def(self*=X());
    differential_vector_class.def("value", &DV::value);
    differential_vector_class.def("jacobian", &DV::jacobian);
    differential_vector_class.def(self_ns::str(self));

    def("compose",(D(*)(const D&,const DV&))&compose);
    def("compose",(DV(*)(const DV&,const DV&))&compose);
}

template void export_differential< Differential<Float> >(const char*);
template void export_differential< Differential<Interval> >(const char*);
template void export_differential< Differential<TaylorModel> >(const char*);

template void export_differential_vector< Differential<Float> >(const char*);
template void export_differential_vector< Differential<Interval> >(const char*);
template void export_differential_vector< Differential<TaylorModel> >(const char*);

void differentiation_submodule() 
{
    //export_differential< DenseDifferential<Float> >("DenseDifferential");
    //export_differential< DenseDifferential<Interval> >("IDenseDifferential");
    //export_differential< SparseDifferential<Float> >("SparseDifferential");
    //export_differential< SparseDifferential<Interval> >("ISparseDifferential");

    //export_differential_vector< DenseDifferential<Float> >("DenseDifferentialVector");
    //export_differential_vector< DenseDifferential<Interval> >("IDenseDifferentialVector");
    //export_differential_vector< SparseDifferential<Float> >("SparseDifferentialVector");
    //export_differential_vector< SparseDifferential<Interval> >("ISpareDifferentialVector");

    export_differential< Differential<Float> >("Differential");
    export_differential< Differential<Interval> >("IntervalDifferential");
    export_differential< Differential<TaylorModel> >("TaylorModelDifferential");

    export_differential_vector< Differential<Float> >("DifferentialVector");
    export_differential_vector< Differential<Interval> >("IntervalDifferentialVector");
    export_differential_vector< Differential<TaylorModel> >("TaylorModelDifferentialVector");
}

