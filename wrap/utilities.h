/***************************************************************************
 *            utilities.h
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

/*! \file utilities.h
 *  Commonly used inline methods for the Python interface.
 */

#ifndef ARIADNE_PYTHON_UTILITIES_H
#define ARIADNE_PYTHON_UTILITIES_H

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>

#include <string>
#include <vector>
#include <set>
#include <map>
#include <sstream>
#include <iostream>

#include "config.h"

#include "utility/array.h"
#include "utility/tuple.h"
#include "utility/container.h"
#include "numeric/numeric.h"
#include "algebra/vector.h"

namespace Ariadne {

class Real;

template<class X> inline const char* python_name(const char* name);
template<> inline const char* python_name<Float>(const char* name) {
    return (StringType("Float")+name).c_str(); }
template<> inline const char* python_name<Integer>(const char* name) {
    return (StringType("Integer")+name).c_str(); }
template<> inline const char* python_name<Rational>(const char* name) {
    return (StringType("Rational")+name).c_str(); }
template<> inline const char* python_name<Real>(const char* name) {
    return (StringType("Real")+name).c_str(); }

template<> inline const char* python_name<ExactFloat>(const char* name) {
    return (StringType("ExactFloat")+name).c_str(); }
template<> inline const char* python_name<ValidatedFloat>(const char* name) {
    return (StringType("ValidatedFloat")+name).c_str(); }
template<> inline const char* python_name<UpperFloat>(const char* name) {
    return (StringType("UpperFloat")+name).c_str(); }
template<> inline const char* python_name<ApproximateFloat>(const char* name) {
    return (StringType("ApproximateFloat")+name).c_str(); }

template<> inline const char* python_name<ExactInterval>(const char* name) {
    return (StringType("ExactInterval")+name).c_str(); }

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


#ifdef HAVE_GMPXX_H
#endif



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
struct to_python< Ariadne::Pair<T1,T2> > {
    to_python() { boost::python::to_python_converter< Ariadne::Pair<T1,T2>, to_python< Ariadne::Pair<T1,T2> > >(); }
    static PyObject* convert(const Ariadne::Pair<T1,T2>& tup) {
        boost::python::list lst;
        lst.append(boost::python::object(tup.first));
        lst.append(boost::python::object(tup.second));
        boost::python::tuple result(lst);
        return boost::python::incref(boost::python::tuple(result).ptr());
    }
    static const PyTypeObject* get_pytype() { return &PyTuple_Type; }
};

template<class T1, class T2>
struct to_python< Ariadne::Tuple<T1,T2> > {
    to_python() { boost::python::to_python_converter< Ariadne::Tuple<T1,T2>, to_python< Ariadne::Tuple<T1,T2> > >(); }
    static PyObject* convert(const Ariadne::Tuple<T1,T2>& tup) {
        boost::python::list lst;
        lst.append(boost::python::object(tup.first));
        lst.append(boost::python::object(tup.second));
        boost::python::tuple result(lst);
        return boost::python::incref(boost::python::tuple(result).ptr());
    }
    static const PyTypeObject* get_pytype() { return &PyTuple_Type; }
};

template<class T1, class T2, class T3>
struct to_python< Ariadne::Tuple<T1,T2,T3> > {
    to_python() { boost::python::to_python_converter< Ariadne::Tuple<T1,T2,T3>, to_python< Ariadne::Tuple<T1,T2,T3> > >(); }
    static PyObject* convert(const Ariadne::Tuple<T1,T2,T3>& tup) {
        boost::python::list lst;
        lst.append(boost::python::object(tup.first));
        lst.append(boost::python::object(tup.second));
        lst.append(boost::python::object(tup.third));
        boost::python::tuple result(lst);
        return boost::python::incref(boost::python::tuple(result).ptr());
    }
    static const PyTypeObject* get_pytype() { return &PyTuple_Type; }
};

template<class T1, class T2, class T3, class T4>
struct to_python< Ariadne::Tuple<T1,T2,T3,T4> > {
    to_python() { boost::python::to_python_converter< Ariadne::Tuple<T1,T2,T3,T4>, to_python< Ariadne::Tuple<T1,T2,T3,T4> > >(); }
    static PyObject* convert(const Ariadne::Tuple<T1,T2,T3,T4>& tup) {
        boost::python::list lst;
        lst.append(boost::python::object(tup.first));
        lst.append(boost::python::object(tup.second));
        lst.append(boost::python::object(tup.third));
        lst.append(boost::python::object(tup.fourth));
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
        for(uint i=0; i!=vec.size(); ++i) {
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
    static void* convertible(PyObject* obj_ptr) {
        if (!PyList_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,boost::python::converter::rvalue_from_python_stage1_data* data) {
        boost::python::list lst = boost::python::extract<boost::python::list>(obj_ptr);
        Array<T> a=Array<T>(len(lst)); for(uint i=0; i!=a.size(); ++i) { (a)[i]=boost::python::extract<T>(lst[i]); }
        void* storage = ((boost::python::converter::rvalue_from_python_storage< Array<T> >*)data)->storage.bytes;
        new (storage) Array<T>(a);
        data->convertible = storage;
    }
};

template<class T>
struct from_python< Ariadne::List<T> > {
    from_python() {
        boost::python::converter::registry::push_back(&convertible,&construct,boost::python::type_id< List<T> >()); }
    static void* convertible(PyObject* obj_ptr) {
        if (!PyList_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,boost::python::converter::rvalue_from_python_stage1_data* data) {
        boost::python::list lst = boost::python::extract<boost::python::list>(obj_ptr);
        List<T> l; l.reserve(len(lst)); for(int i=0; i!=len(lst); ++i) { l.append(boost::python::extract<T>(lst[i])); }
        void* storage = ((boost::python::converter::rvalue_from_python_storage< Array<T> >*)data)->storage.bytes;
        new (storage) List<T>(l);
        data->convertible = storage;
    }
};

template<class T>
struct from_python_list< Ariadne::Array<T> > {
    from_python_list() {
        boost::python::converter::registry::push_back(&convertible,&construct,boost::python::type_id< Array<T> >()); }
    static void* convertible(PyObject* obj_ptr) {
        if (!PyList_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,boost::python::converter::rvalue_from_python_stage1_data* data) {
        boost::python::list lst = boost::python::extract<boost::python::list>(obj_ptr);
        Array<T> a=Array<T>(len(lst)); for(uint i=0; i!=a.size(); ++i) { (a)[i]=boost::python::extract<T>(lst[i]); }
        void* storage = ((boost::python::converter::rvalue_from_python_storage< Array<T> >*)data)->storage.bytes;
        new (storage) Array<T>(a);
        data->convertible = storage;
    }
};

template<class T>
struct from_python_tuple< Ariadne::Array<T> > {
    from_python_tuple() {
        boost::python::converter::registry::push_back(&convertible,&construct,boost::python::type_id< Array<T> >()); }
    static void* convertible(PyObject* obj_ptr) {
        if (!PyTuple_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,boost::python::converter::rvalue_from_python_stage1_data* data) {
        boost::python::tuple tup = boost::python::extract<boost::python::tuple>(obj_ptr);
        Array<T> a(len(tup)); for(uint i=0; i!=a.size(); ++i) { (a)[i]=boost::python::extract<T>(tup[i]); }
        void* storage = ((boost::python::converter::rvalue_from_python_storage< Array<T> >*)data)->storage.bytes;
        new (storage) Array<T>(a);
        data->convertible = storage;
    }
};





template<class C, class I0, class I1, class S> inline
S __getslice__(const C& c, const I0& i0, const I1& i1) { return project(c,range(i0,i1)); }

template<class C, class I, class X> inline
X __getitem__(const C& c, const I& i) { return c[i]; }

template<class C, class I, class X> inline
void __setitem__(C& c, const I& i, const X& x) { c[i]=x; }

template<class C, class I, class J, class X> inline
X __getitem2__(const C& c, const I& i, const J& j) { return c[i][j]; }

template<class C, class I, class J, class X> inline
void __setitem2__(C& c, const I& i, const J& j, const X& x) { c[i][j]=x; }



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
R __and__(const A1& a1, const A2& a2) { return static_cast<R>(a1 && a2); }

template<class R, class A1, class A2>
R __or__(const A1& a1, const A2& a2) { return static_cast<R>(a1 || a2); }

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
void export_array(const char* name)
{
    boost::python::class_< Array<T> > array_class(name,boost::python::init < Array<T> >());
    array_class.def(boost::python::init<int>());
    array_class.def("__len__", &Array<T>::size);
    array_class.def("__getitem__",&__getitem__< Array<T>, uint, T>);
    array_class.def("__setitem__",&__setitem__< Array<T>, uint, T>);
    array_class.def(boost::python::self_ns::str(boost::python::self));
}




}

#endif /* ARIADNE_PYTHON_UTILITIES_H */
