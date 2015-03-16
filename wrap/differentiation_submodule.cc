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

#include "boost_python.h"
#include "utilities.h"

#include <boost/python.hpp>

#include "utility/array.h"
#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "algebra/covector.h"
#include "algebra/differential.h"

using namespace boost::python;
using namespace Ariadne;

template<class X> Void read_array(Array<X>&, const boost::python::object& obj) { }
inline Nat compute_polynomial_data_size(Nat rs, Nat as, Nat d) { return rs*Ariadne::bin(d+as,as); }

namespace Ariadne {

// FIXME: Ensure all valid arithmetic and comparisons are defined!
inline auto operator==(BoundedFloat64 x, Int n) -> decltype(x==ExactFloat64(n)) { return x==ExactFloat64(n); }
inline auto operator!=(BoundedFloat64 x, Int n) -> decltype(x!=ExactFloat64(n)) { return x!=ExactFloat64(n); }
inline auto operator> (BoundedFloat64 x, Int n) -> decltype(x> ExactFloat64(n)) { return x> ExactFloat64(n); }
inline auto operator*=(ApproximateFloat64 x, Int n) -> decltype(x*=ApproximateFloat64(n)) { return x*=ApproximateFloat64(n); }

template<class X>
struct to_python_dict< Ariadne::Expansion<X>  > {
    to_python_dict() { boost::python::to_python_converter< Ariadne::Expansion<X>, to_python_dict< Ariadne::Expansion<X> > >(); }
    static PyObject* convert(const Ariadne::Expansion<X>& e) {
        Nat n=e.argument_size();
        boost::python::dict res;
        boost::python::list lst;
        for(Nat i=0; i!=n; ++i) { lst.append(0); }
        Ariadne::MultiIndex a;
        X c;
        for(typename Expansion<X>::ConstIterator iter=e.begin(); iter!=e.end(); ++iter) {
            a=iter->key();
            c=iter->data();
            for(Nat i=0; i!=a.size(); ++i) { Int ai=a[i]; lst[i]=ai; }
            boost::python::tuple tup(lst);
            //res[tup]=boost::python::object(c);
            res[boost::python::object(a)]=boost::python::object(c);
        }
        return boost::python::incref(boost::python::dict(res).ptr());
    }
    static const PyTypeObject* get_pytype() { return &PyDict_Type; }
};

template<class X>
struct to_python_list< Ariadne::Expansion<X>  > {
    to_python_list() { boost::python::to_python_converter< Ariadne::Expansion<X>, to_python_list< Ariadne::Expansion<X> > >(); }
    static PyObject* convert(const Ariadne::Expansion<X>& e) {
        Nat n=e.argument_size();
        boost::python::list res;
        boost::python::list alst;
        for(Nat i=0; i!=n; ++i) { alst.append(0); }
        std::cerr<<"Here\n";
        boost::python::list pr; pr.append(0); pr.append(0);
        Ariadne::MultiIndex a;
        X c;
        for(typename Expansion<X>::ConstIterator iter=e.begin(); iter!=e.end(); ++iter) {
            a=iter->key();
            c=iter->data();
            for(Nat i=0; i!=n; ++i) { Int ai=a[i]; alst[i]=ai; }
            pr[0]=boost::python::tuple(alst);
            pr[1]=boost::python::object(c);
            res.append(boost::python::tuple(pr));
        }
        return boost::python::incref(boost::python::list(res).ptr());
    }
    static const PyTypeObject* get_pytype() { return &PyList_Type; }
};

/*
   static PyObject* convert(const Tuple<T1,T2,T3,T4>& tup) {
        boost::python::list lst;
        lst.append(boost::python::object(tup.first));
        lst.append(boost::python::object(tup.second));
        lst.append(boost::python::object(tup.third));
        lst.append(boost::python::object(tup.fourth));
        boost::python::tuple result(lst);
        return boost::python::incref(boost::python::tuple(result).ptr());
    static PyObject* convert(const std::map<K,V>& map) {
        boost::python::dict result;
        for(typename std::map<K,V>::ConstIterator iter=map.begin(); iter!=map.end(); ++iter) {
            result[boost::python::object(iter->first)]=boost::python::object(iter->second);
        }
        return boost::python::incref(boost::python::dict(result).ptr());
*/
}


template<class DIFF>
DIFF*
make_differential(const Nat& as, const Nat& d, const boost::python::object& obj)
{
    typedef typename DIFF::ValueType X;
    DIFF* result=new DIFF(as,d);
    Array<X> data;
    read_array(data,obj);
    std::cerr<<"polynomial_data_size("<<as<<","<<d<<")="<<compute_polynomial_data_size(1u,as,d)<<"\n";
    std::cerr<<"data="<<data<<"\n";
    assert(data.size()==compute_polynomial_data_size(1u,as,d));
    MultiIndex i(as);
    const X* ptr=data.begin();
    while(i.degree()<=d) {
        //result[i]=*ptr; ++i; ++ptr;
        result->expansion().append(i,*ptr); ++i; ++ptr;
    }
    return result;
}

template<class DIFF>
DIFF*
make_sparse_differential(const boost::python::object& obj,const Nat& d)
{
    typedef typename DIFF::ValueType X;
    Expansion<X> expansion = boost::python::extract< Expansion<X> >(obj);
    DIFF* result=new DIFF(expansion,d);
    return result;
}


template<class DIFF>
boost::python::list
make_differential_variables(const Nat& d, const Vector<typename DIFF::NumericType>& x)
{
    boost::python::list result;
    for(Nat i=0; i!=x.size(); ++i) {
        result.append(DIFF::variable(x.size(),d,numeric_cast<typename DIFF::NumericType>(x[i]),i));
    }
    return result;
}


template<class DIFF>
Vector<DIFF>*
make_differential_vector(const Nat& rs, const Nat& as, const Nat& d, const boost::python::object& obj)
{
    typedef typename DIFF::ValueType X;
    Array<X> data;
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
Void set_item(C& c, const I& i, const X& x) { c[i]=x; }

template<class C, class I, class J, class X> inline
Void matrix_set_item(C& c, const I& i, const J& j, const X& x) { c[i][j]=x; }


namespace Ariadne {

template<class X> OutputStream& operator<<(OutputStream& os, const PythonRepresentation< Expansion<X> >& repr);

template<class X> OutputStream& operator<<(OutputStream& os, const PythonRepresentation< Differential<X> >& repr) {
    const Differential<X>& diff=repr.reference();
    os << python_name<X>("Differential") << "(" << python_representation(diff.expansion()) << "," << diff.degree() << ")";
    //os << python_name<X>("Differential") << "(" << diff.argument_size() << "," << diff.degree() << "," << python_representation(diff.expansion()) << ")";
    return os;
}

}



template<class DIFF>
Void export_differential(const char* name)
{
    typedef typename DIFF::ValueType X;
    typedef Vector<X> V;
    typedef Series<X> S;
    typedef DIFF D;
    typedef Vector<D> DV;


    class_<D> differential_class(name, init<D>() );
    //differential_class.def("__init__", make_constructor(&make_differential<D>) );
    differential_class.def("__init__", make_constructor(&make_sparse_differential<D>) );
    //differential_class.def( init< D >());
    differential_class.def( init< Nat, Nat >());
    differential_class.def( init< Expansion<X>, Nat >());
    differential_class.def("__getitem__", &get_item<D,MultiIndex,X>);
    differential_class.def("__setitem__",&set_item<D,MultiIndex,X>);
    differential_class.def(-self);
    differential_class.def(self+self);
    differential_class.def(self-self);
    differential_class.def(self*self);
    differential_class.def(self/self);
    differential_class.def(self+X());
    differential_class.def(self-X());
    differential_class.def(self*X());
    differential_class.def(self/X());
    differential_class.def(X()+self);
    differential_class.def(X()-self);
    differential_class.def(X()*self);
    differential_class.def(self+=self);
    differential_class.def(self-=self);
    //differential_class.def(self*=self);
    differential_class.def(self+=X());
    differential_class.def(self-=X());
    differential_class.def(self*=X());
    differential_class.def(self/=X());
    differential_class.def(self_ns::str(self));
    //differential_class.def("__repr__", (StringType(*)(const D&)) &__repr__);
    differential_class.def("__repr__", &__repr__<D>);

    differential_class.def("value",&D::value,return_value_policy<copy_const_reference>());
    differential_class.def("gradient",(Covector<X>(D::*)()const)&D::gradient);
    differential_class.def("hessian", (Matrix<X>(D::*)()const)&D::hessian);
    differential_class.def("expansion", (Expansion<X>const&(D::*)()const)&D::expansion,return_value_policy<copy_const_reference>());

    differential_class.def("constant",(D(*)(SizeType, DegreeType, const X&))&D::constant);
    differential_class.def("variable",(D(*)(SizeType, DegreeType, const X&, SizeType))&D::variable);
    differential_class.def("variables",(Vector<D>(*)(DegreeType, const Vector<X>&))&D::variables);

    differential_class.staticmethod("constant");
    differential_class.staticmethod("variable");
    differential_class.staticmethod("variables");

    def("derivative", (D(*)(const D&, SizeType))&derivative);
    def("antiderivative", (D(*)(const D&, SizeType))&antiderivative);

    def("pos",(D(*)(const D&))&pos);
    def("neg",(D(*)(const D&))&neg);
    def("rec",(D(*)(const D&))&rec);
    def("pow",(D(*)(const D&, Int))&pow);

    def("sqrt", (D(*)(const D&))&sqrt);
    def("exp", (D(*)(const D&))&exp);
    def("log", (D(*)(const D&))&log);
    def("sin", (D(*)(const D&))&sin);
    def("cos", (D(*)(const D&))&cos);
    def("tan", (D(*)(const D&))&tan);
    def("atan", (D(*)(const D&))&atan);
}

template<class DIFF>
Void
export_differential_vector(const char* name)
{
    typedef typename DIFF::ValueType X;
    typedef Vector<X> V;
    typedef Series<X> S;
    typedef DIFF D;
    typedef Vector<D> DV;

    class_<DV> differential_vector_class(name, init<DV>());
    differential_vector_class.def("__init__", make_constructor(&make_differential_vector<D>) );
    differential_vector_class.def( init< Nat, Nat, Nat >());
    differential_vector_class.def("__getitem__", &matrix_get_item<DV,Int,MultiIndex,X>);
    differential_vector_class.def("__getitem__", &get_item<DV,Int,D>);
    differential_vector_class.def("__setitem__",&set_item<DV,Int,X>);
    differential_vector_class.def("__setitem__",&set_item<DV,Int,D>);
    differential_vector_class.def("__neg__",&__neg__<DV,DV>);
    differential_vector_class.def("__add__",&__add__<DV,DV,DV>);
    differential_vector_class.def("__sub__",&__sub__<DV,DV,DV>);
    differential_vector_class.def("__add__",&__add__<DV,DV,V>);
    differential_vector_class.def("__sub__",&__sub__<DV,DV,V>);
    //differential_vector_class.def("__mul__",&__mul__<DV,DV,D>);
    //differential_vector_class.def("__div__",&__div__<DV,DV,D>);
    differential_vector_class.def("__rmul__",&__rmul__<DV,DV,X>);
    differential_vector_class.def("__mul__",&__mul__<DV,DV,X>);
    differential_vector_class.def("__div__",&__div__<DV,DV,X>);
    differential_vector_class.def("value", &DV::value);
    differential_vector_class.def("jacobian", &DV::jacobian);
    differential_vector_class.def(self_ns::str(self));

    def("compose",(D(*)(const D&,const DV&))&compose);
    def("compose",(DV(*)(const DV&,const DV&))&compose);

    def("lie_derivative", (DV(*)(const DV&,const DV&))&lie_derivative);
}

template Void export_differential< Differential<ApproximateFloat64> >(const char*);
template Void export_differential< Differential<BoundedFloat64> >(const char*);

template Void export_differential_vector< Differential<ApproximateFloat64> >(const char*);
template Void export_differential_vector< Differential<BoundedFloat64> >(const char*);

Void differentiation_submodule()
{
    to_python_dict < Expansion<ApproximateFloat64> >();
    to_python_dict < Expansion<BoundedFloat64> >();

    export_differential< Differential<ApproximateFloat64> >("ApproximateDifferential");
    export_differential< Differential<BoundedFloat64> >("ValidatedDifferential");

    export_differential_vector< Differential<ApproximateFloat64> >("ApproximateDifferentialVector");
    export_differential_vector< Differential<BoundedFloat64> >("ValidatedDifferentialVector");
}

