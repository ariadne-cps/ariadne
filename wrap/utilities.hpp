/***************************************************************************
 *            utilities.hpp
 *
 *  Copyright  2005-8  Alberto Casagrande, Pieter Collins
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

/*! \file utilities.hpp
 *  Commonly used inline methods for the Python interface.
 */

#ifndef ARIADNE_PYTHON_UTILITIES_HPP
#define ARIADNE_PYTHON_UTILITIES_HPP

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>

#include <string>
#include <vector>
#include <set>
#include <map>
#include <sstream>
#include <iostream>

#include "config.h"

#include "utility/array.hpp"
#include "utility/tuple.hpp"
#include "utility/container.hpp"
#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"

namespace Ariadne {

class Real;

template<class X> inline std::string python_name(const std::string& name) {
    return class_name<X>()+name; }

template<class T> struct to_python;
template<class T> struct to_python_list;
template<class T> struct to_python_dict;
template<class T> struct to_python_tuple;
template<class T> struct to_python_set;

template<class T> struct from_python;
template<class T> struct from_python_list;
template<class T> struct from_python_tuple;
template<class T> struct from_python_dict;
template<class T> struct from_python_set;




template<class T1, class T2>
struct to_python< std::pair<T1,T2> > {
    to_python() { boost::python::to_python_converter< std::pair<T1,T2>, to_python< std::pair<T1,T2> > >(); }
    static PyObject* convert(const std::pair<T1,T2>& pair) {
        boost::python::list lst;
        lst.append(boost::python::object(pair.first));
        lst.append(boost::python::object(pair.second));
        boost::python::tuple result(lst);
        return boost::python::incref(boost::python::tuple(result).ptr());
    }
    static const PyTypeObject* get_pytype() { return &PyTuple_Type; }
};


template<class T>
struct to_python< std::vector<T> > {
    to_python() { boost::python::to_python_converter< std::vector<T>, to_python< std::vector<T> > >(); }

    static PyObject* convert(const std::vector<T>& vec) {
        boost::python::list result;
        for(typename std::vector<T>::ConstIterator iter=vec.begin(); iter!=vec.end(); ++iter) {
            result.append(boost::python::object(*iter));
        }
        return boost::python::incref(boost::python::list(result).ptr());
    }
    static const PyTypeObject* get_pytype() { return &PyList_Type; }
};


template<class T>
struct to_python<std::set<T> > {
    to_python() { boost::python::to_python_converter< std::set<T>, to_python< std::set<T> > >(); }
    static PyObject* convert(const std::set<T>& set) {
        boost::python::list values;
        for(typename std::set<T>::ConstIterator iter=set.begin(); iter!=set.end(); ++iter) {
            values.append(boost::python::object(*iter));
        }
        PyObject* result=PySet_New(values.ptr());
        return result;
    }
    static const PyTypeObject* get_pytype() { return &PySet_Type; }
};

template<class T>
struct to_python_list<std::set<T> > {
    to_python_list() { boost::python::to_python_converter< std::set<T>, to_python_list< std::set<T> > >(); }
    static PyObject* convert(const std::set<T>& set) {
        boost::python::list result;
        for(typename std::set<T>::ConstIterator iter=set.begin(); iter!=set.end(); ++iter) {
            result.append(boost::python::object(*iter));
        }
        return boost::python::incref(boost::python::list(result).ptr());
    }
    static const PyTypeObject* get_pytype() { return &PyList_Type; }
};

template<class K, class V>
struct to_python<std::map<K,V> > {
    to_python() { boost::python::to_python_converter< std::map<K,V>, to_python< std::map<K,V> > >(); }
    static PyObject* convert(const std::map<K,V>& map) {
        boost::python::dict result;
        for(typename std::map<K,V>::ConstIterator iter=map.begin(); iter!=map.end(); ++iter) {
            result[boost::python::object(iter->first)]=boost::python::object(iter->second);
        }
        return boost::python::incref(boost::python::dict(result).ptr());
    }
    static const PyTypeObject* get_pytype() { return &PyDict_Type; }
};


template<class T1, class T2>
struct to_python< std::tuple<T1,T2> > {
    to_python() { boost::python::to_python_converter< Ariadne::Tuple<T1,T2>, to_python< Ariadne::Tuple<T1,T2> > >(); }
    static PyObject* convert(const Ariadne::Tuple<T1,T2>& tup) {
        boost::python::list lst;
        lst.append(boost::python::object(std::get<0>(tup)));
        lst.append(boost::python::object(std::get<1>(tup)));
        boost::python::tuple result(lst);
        return boost::python::incref(boost::python::tuple(result).ptr());
    }
    static const PyTypeObject* get_pytype() { return &PyTuple_Type; }
};

template<class T1, class T2, class T3>
struct to_python< std::tuple<T1,T2,T3> > {
    to_python() { boost::python::to_python_converter< Ariadne::Tuple<T1,T2,T3>, to_python< Ariadne::Tuple<T1,T2,T3> > >(); }
    static PyObject* convert(const Ariadne::Tuple<T1,T2,T3>& tup) {
        boost::python::list lst;
        lst.append(boost::python::object(std::get<0>(tup)));
        lst.append(boost::python::object(std::get<1>(tup)));
        lst.append(boost::python::object(std::get<2>(tup)));
        boost::python::tuple result(lst);
        return boost::python::incref(boost::python::tuple(result).ptr());
    }
    static const PyTypeObject* get_pytype() { return &PyTuple_Type; }
};

template<class T1, class T2, class T3, class T4>
struct to_python< std::tuple<T1,T2,T3,T4> > {
    to_python() { boost::python::to_python_converter< Ariadne::Tuple<T1,T2,T3,T4>, to_python< Ariadne::Tuple<T1,T2,T3,T4> > >(); }
    static PyObject* convert(const Ariadne::Tuple<T1,T2,T3,T4>& tup) {
        boost::python::list lst;
        lst.append(boost::python::object(std::get<0>(tup)));
        lst.append(boost::python::object(std::get<1>(tup)));
        lst.append(boost::python::object(std::get<2>(tup)));
        lst.append(boost::python::object(std::get<3>(tup)));
        boost::python::tuple result(lst);
        return boost::python::incref(boost::python::tuple(result).ptr());
    }
    static const PyTypeObject* get_pytype() { return &PyTuple_Type; }
};

template<class T>
struct to_python< Ariadne::Array<T> > {
    to_python() { boost::python::to_python_converter< Ariadne::Array<T>, to_python< Ariadne::Array<T> > >(); }
    static PyObject* convert(const Ariadne::Array<T>& ary) {
        boost::python::list result;
        for(typename Ariadne::Array<T>::ConstIterator iter=ary.begin(); iter!=ary.end(); ++iter) {
            result.append(boost::python::object(*iter));
        }
        return boost::python::incref(boost::python::list(result).ptr());
    }
    static const PyTypeObject* get_pytype() { return &PyList_Type; }
};

template<class T>
struct to_python< Ariadne::List<T> > {
    to_python() { boost::python::to_python_converter< Ariadne::List<T>, to_python< Ariadne::List<T> > >(); }
    static PyObject* convert(const Ariadne::List<T>& lst) {
        boost::python::list result;
        for(typename List<T>::ConstIterator iter=lst.begin(); iter!=lst.end(); ++iter) {
            result.append(boost::python::object(*iter));
        }
        return boost::python::incref(boost::python::list(result).ptr());
    }
    static const PyTypeObject* get_pytype() { return &PyList_Type; }
};

template<class T>
struct to_python< Ariadne::Set<T> > {
    to_python() { boost::python::to_python_converter< Ariadne::Set<T>, to_python< Ariadne::Set<T> > >(); }
    static PyObject* convert(const Ariadne::Set<T>& set) {
        boost::python::list values;
        for(typename Ariadne::Set<T>::ConstIterator iter=set.begin(); iter!=set.end(); ++iter) {
            values.append(boost::python::object(*iter));
        }
        PyObject* result=PySet_New(values.ptr());
        return result;
    }
    static const PyTypeObject* get_pytype() { return &PySet_Type; }
};

template<class T>
struct to_python_list< Ariadne::Set<T> > {
    to_python_list() { boost::python::to_python_converter< Ariadne::Set<T>, to_python_list< Ariadne::Set<T> > >(); }
    static PyObject* convert(const Set<T>& set) {
        boost::python::list result;
        for(typename Set<T>::ConstIterator iter=set.begin(); iter!=set.end(); ++iter) {
            result.append(boost::python::object(*iter));
        }
        return boost::python::incref(boost::python::list(result).ptr());
    }
    static const PyTypeObject* get_pytype() { return &PyList_Type; }
};

template<class K, class V>
struct to_python<Ariadne::Map<K,V> > {
    to_python() { boost::python::to_python_converter< Ariadne::Map<K,V>, to_python< Ariadne::Map<K,V> > >(); }
    static PyObject* convert(const Ariadne::Map<K,V>& map) {
        boost::python::dict result;
        for(typename Ariadne::Map<K,V>::ConstIterator iter=map.begin(); iter!=map.end(); ++iter) {
            result[boost::python::object(iter->first)]=boost::python::object(iter->second);
        }
        return boost::python::incref(boost::python::dict(result).ptr());
    }
    static const PyTypeObject* get_pytype() { return &PyDict_Type; }
};


template<class X>
struct to_python< Vector<X> > {
    to_python() { boost::python::to_python_converter< Vector<X>, to_python< Vector<X> > >(); }
    static PyObject* convert(const Vector<X>& vec) {
        boost::python::list result;
        for(Nat i=0; i!=vec.size(); ++i) {
            result.append(boost::python::object(vec[i]));
        }
        return boost::python::incref(boost::python::list(result).ptr());
    }
    static const PyTypeObject* get_pytype() { return &PyList_Type; }
};


template<class T>
struct from_python< Ariadne::Array<T> > {
    from_python() {
        boost::python::converter::registry::push_back(&convertible,&construct,boost::python::type_id< Array<T> >()); }
    static Void* convertible(PyObject* obj_ptr) {
        if (!PyList_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,boost::python::converter::rvalue_from_python_stage1_data* data) {
        boost::python::list lst = boost::python::extract<boost::python::list>(obj_ptr);
        Array<T> a=Array<T>(len(lst)); for(Nat i=0; i!=a.size(); ++i) { (a)[i]=boost::python::extract<T>(lst[i]); }
        Void* storage = ((boost::python::converter::rvalue_from_python_storage< Array<T> >*)data)->storage.bytes;
        new (storage) Array<T>(a);
        data->convertible = storage;
    }
};

template<class T>
struct from_python< Ariadne::List<T> > {
    from_python() {
        boost::python::converter::registry::push_back(&convertible,&construct,boost::python::type_id< List<T> >()); }
    static Void* convertible(PyObject* obj_ptr) {
        if (!PyList_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,boost::python::converter::rvalue_from_python_stage1_data* data) {
        boost::python::list lst = boost::python::extract<boost::python::list>(obj_ptr);
        List<T> l; l.reserve(len(lst)); for(Int i=0; i!=len(lst); ++i) { l.append(boost::python::extract<T>(lst[i])); }
        Void* storage = ((boost::python::converter::rvalue_from_python_storage< Array<T> >*)data)->storage.bytes;
        new (storage) List<T>(l);
        data->convertible = storage;
    }
};

template<class T>
struct from_python_list< Ariadne::Array<T> > {
    from_python_list() {
        boost::python::converter::registry::push_back(&convertible,&construct,boost::python::type_id< Array<T> >()); }
    static Void* convertible(PyObject* obj_ptr) {
        if (!PyList_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,boost::python::converter::rvalue_from_python_stage1_data* data) {
        boost::python::list lst = boost::python::extract<boost::python::list>(obj_ptr);
        Array<T> a=Array<T>(len(lst)); for(Nat i=0; i!=a.size(); ++i) { (a)[i]=boost::python::extract<T>(lst[i]); }
        Void* storage = ((boost::python::converter::rvalue_from_python_storage< Array<T> >*)data)->storage.bytes;
        new (storage) Array<T>(a);
        data->convertible = storage;
    }
};

template<class T>
struct from_python_tuple< Ariadne::Array<T> > {
    from_python_tuple() {
        boost::python::converter::registry::push_back(&convertible,&construct,boost::python::type_id< Array<T> >()); }
    static Void* convertible(PyObject* obj_ptr) {
        if (!PyTuple_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,boost::python::converter::rvalue_from_python_stage1_data* data) {
        boost::python::tuple tup = boost::python::extract<boost::python::tuple>(obj_ptr);
        Array<T> a(len(tup)); for(Nat i=0; i!=a.size(); ++i) { (a)[i]=boost::python::extract<T>(tup[i]); }
        Void* storage = ((boost::python::converter::rvalue_from_python_storage< Array<T> >*)data)->storage.bytes;
        new (storage) Array<T>(a);
        data->convertible = storage;
    }
};





template<class C, class I0, class I1, class S> inline
S __getslice__(const C& c, const I0& i0, const I1& i1) { return project(c,range(i0,i1)); }

template<class C, class I, class X> inline
X __getitem__(const C& c, const I& i) { return c[i]; }

template<class C, class I, class X> inline
Void __setitem__(C& c, const I& i, const X& x) { c[i]=x; }

template<class C, class I, class J, class X> inline
X __getitem2__(const C& c, const I& i, const J& j) { return c[i][j]; }

template<class C, class I, class J, class X> inline
Void __setitem2__(C& c, const I& i, const J& j, const X& x) { c[i][j]=x; }



template<class R, class A>
R __pos__(const A& a) { return static_cast<R>(a); }

template<class R, class A>
R __neg__(const A& a) { return static_cast<R>(-a); }

template<class R, class A1, class A2>
R __add__(const A1& a1, const A2& a2) { return static_cast<R>(a1+a2); }

template<class R, class A1, class A2>
R __sub__(const A1& a1, const A2& a2) { return static_cast<R>(a1-a2); }

template<class R, class A1, class A2>
R __mul__(const A1& a1, const A2& a2) { return static_cast<R>(a1*a2); }

template<class R, class A1, class A2>
R __div__(const A1& a1, const A2& a2) { return static_cast<R>(a1/a2); }

template<class R, class A1, class A2>
R __radd__(const A1& a1, const A2& a2) { return static_cast<R>(a2+a1); }

template<class R, class A1, class A2>
R __rsub__(const A1& a1, const A2& a2) { return static_cast<R>(a2-a1); }

template<class R, class A1, class A2>
R __rmul__(const A1& a1, const A2& a2) { return static_cast<R>(a2*a1); }

template<class R, class A1, class A2>
R __rdiv__(const A1& a1, const A2& a2) { return static_cast<R>(a2/a1); }

template<class R, class A1, class A2>
R __pow__(const A1& a1, const A2& a2) { return static_cast<R>(pow(a1,a2)); }

template<class R, class A1, class A2>
R __and__(const A1& a1, const A2& a2) { return static_cast<R>(a1 && a2); }

template<class R, class A1, class A2>
R __or__(const A1& a1, const A2& a2) { return static_cast<R>(a1 || a2); }

template<class R, class A1, class A2>
R __bitor__(const A1& a1, const A2& a2) { return static_cast<R>(a1 | a2); }

template<class R, class A>
R __not__(const A& a) { return static_cast<R>(!a); }

template<class R, class A1, class A2>
R __eq__(const A1& a1, const A2& a2) { return static_cast<R>(a1==a2); }

template<class R, class A1, class A2>
R __ne__(const A1& a1, const A2& a2) { return static_cast<R>(a1!=a2); }

template<class R, class A1, class A2>
R __gt__(const A1& a1, const A2& a2) { return static_cast<R>(a1>a2); }

template<class R, class A1, class A2>
R __lt__(const A1& a1, const A2& a2) { return static_cast<R>(a1<a2); }

template<class R, class A1, class A2>
R __ge__(const A1& a1, const A2& a2) { return static_cast<R>(a1>=a2); }

template<class R, class A1, class A2>
R __le__(const A1& a1, const A2& a2) { return static_cast<R>(a1<=a2); }


template<class A> auto _pos_(const A& a) -> decltype(pos(a)) { return pos(a); }
template<class A> auto _neg_(const A& a) -> decltype(neg(a)) { return neg(a); }
template<class A> auto _rec_(const A& a) -> decltype(rec(a)) { return rec(a); }
template<class A> auto _sqr_(const A& a) -> decltype(sqr(a)) { return sqr(a); }
template<class A, class N> auto _pow_(const A& a, N n) -> decltype(pow(a,n)) { return pow(a,n); }
template<class A> auto _sqrt_(const A& a) -> decltype(sqrt(a)) { return sqrt(a); }
template<class A> auto _exp_(const A& a) -> decltype(exp(a)) { return exp(a); }
template<class A> auto _log_(const A& a) -> decltype(log(a)) { return log(a); }
template<class A> auto _sin_(const A& a) -> decltype(sin(a)) { return sin(a); }
template<class A> auto _cos_(const A& a) -> decltype(cos(a)) { return cos(a); }
template<class A> auto _tan_(const A& a) -> decltype(tan(a)) { return tan(a); }
template<class A> auto _atan_(const A& a) -> decltype(atan(a)) { return atan(a); }
template<class A1, class A2> auto _max_(const A1& a1, const A2& a2) -> decltype(max(a1,a2)) { return max(a1,a2); }
template<class A1, class A2> auto _min_(const A1& a1, const A2& a2) -> decltype(min(a1,a2)) { return min(a1,a2); }
template<class A> auto _abs_(const A& a) -> decltype(abs(a)) { return abs(a); }

template<class... AS> inline auto _evaluate_(AS... as) -> decltype(evaluate(as...)) { return evaluate(as...); }
template<class... AS> inline auto _partial_evaluate_(AS... as) -> decltype(partial_evaluate(as...)) { return partial_evaluate(as...); }
template<class... AS> inline auto _unchecked_evaluate_(AS... as) -> decltype(unchecked_evaluate(as...)) { return unchecked_evaluate(as...); }
template<class... AS> inline auto _compose_(AS... as) -> decltype(compose(as...)) { return compose(as...); }
template<class... AS> inline auto _unchecked_compose_(AS... as) -> decltype(unchecked_compose(as...)) { return unchecked_compose(as...); }

template<class... AS> inline auto _join_(AS... as) -> decltype(join(as...)) { return join(as...); }
template<class... AS> inline auto _combine_(AS... as) -> decltype(combine(as...)) { return combine(as...); }

template<class... AS> inline auto _midpoint_(AS... as) -> decltype(midpoint(as...)) { return midpoint(as...); }
template<class... AS> inline auto _embed_(AS... as) -> decltype(embed(as...)) { return embed(as...); }
template<class... AS> inline auto _extension_(AS... as) -> decltype(extension(as...)) { return extension(as...); }
template<class... AS> inline auto _restriction_(AS... as) -> decltype(restriction(as...)) { return restriction(as...); }
    template<class... AS> inline auto _split_(AS... as) -> decltype(split(as...)) { return split(as...); }
template<class... AS> inline auto _derivative_(AS... as) -> decltype(derivative(as...)) { return derivative(as...); }
template<class... AS> inline auto _antiderivative_(AS... as) -> decltype(antiderivative(as...)) { return antiderivative(as...); }

template<class... AS> inline auto _refinement_(AS... as) -> decltype(refinement(as...)) { return refinement(as...); }
template<class... AS> inline auto _refines_(AS... as) -> decltype(refines(as...)) { return refines(as...); }
template<class... AS> inline auto _inconsistent_(AS... as) -> decltype(inconsistent(as...)) { return inconsistent(as...); }

template<class T> StringType __cstr__(const T& t) {
    StringStream ss; ss << t; return ss.str(); }

template<class T> StringType __crepr__(const T& t) {
    StringStream ss; ss << representation(t); return ss.str(); }


template<class T> struct PythonRepresentation {
    const T* pointer;
    PythonRepresentation(const T& t) : pointer(&t) { }
    const T& reference() const { return *pointer; }
};

template<class T> PythonRepresentation<T>
python_representation(const T& t) {
    return PythonRepresentation<T>(t); }

template<class T> StringType __repr__(const T& t) {
    StringStream ss;
    ss << PythonRepresentation<T>(t);
    return ss.str();
}

template<class T>
Void export_array(const char* name)
{
    boost::python::class_< Array<T> > array_class(name,boost::python::init < Array<T> >());
    array_class.def(boost::python::init<Int>());
    array_class.def("__len__", &Array<T>::size);
    array_class.def("__getitem__",&__getitem__< Array<T>, Nat, T>);
    array_class.def("__setitem__",&__setitem__< Array<T>, Nat, T>);
    array_class.def(boost::python::self_ns::str(boost::python::self));
}




} // namespace Ariadne

#endif /* ARIADNE_PYTHON_UTILITIES_HPP */
