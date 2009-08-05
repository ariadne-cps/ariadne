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

#include "array.h"
#include "tuple.h"
#include "container.h"
#include "numeric.h"

namespace Ariadne {

template<class X> inline const char* python_name(const char* name);
template<> inline const char* python_name<Float>(const char* name) {
    return (std::string("")+name).c_str(); }
template<> inline const char* python_name<Interval>(const char* name) {
    return (std::string("I")+name).c_str(); }


template<class T> std::string __cstr__(const T& t) {
    std::stringstream ss;
    ss << t;
    return ss.str();
}



void read(bool&, const boost::python::object&);
void read(int&, const boost::python::object&);
void read(long int&, const boost::python::object&);
void read(unsigned int&, const boost::python::object&);
void read(unsigned long int&, const boost::python::object&);
void read(double&, const boost::python::object&);
void read(Float&, const boost::python::object&);
void read(Interval&, const boost::python::object&);

#ifdef HAVE_GMPXX_H
//void read(Integer&, const boost::python::object&);
void read(Rational&, const boost::python::object&);
#endif


template<class T1, class T2>
struct std_pair_to_python {
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
struct std_vector_to_python {
    static PyObject* convert(const std::vector<T>& lst) {
        boost::python::list result;
        for(typename std::vector<T>::const_iterator iter=lst.begin(); iter!=lst.end(); ++iter) {
            result.append(boost::python::object(*iter));
        }
        return boost::python::incref(boost::python::list(result).ptr());
    }
    static const PyTypeObject* get_pytype() { return &PyList_Type; }
};

template<class T>
struct std_set_to_python {
    static PyObject* convert(const std::set<T>& set) {
        boost::python::list values;
        for(typename std::set<T>::const_iterator iter=set.begin(); iter!=set.end(); ++iter) {
            values.append(boost::python::object(*iter));
        }
        PyObject* result=PySet_New(values.ptr());
        return result;
    }
    static const PyTypeObject* get_pytype() { return &PySet_Type; }
};

template<class T>
struct std_set_to_python_list {
    static PyObject* convert(const std::set<T>& set) {
        boost::python::list result;
        for(typename std::set<T>::const_iterator iter=set.begin(); iter!=set.end(); ++iter) {
            result.append(boost::python::object(*iter));
        }
        return boost::python::incref(boost::python::list(result).ptr());
    }
    static const PyTypeObject* get_pytype() { return &PyList_Type; }
};

template<class K, class V>
struct std_map_to_python {
    static PyObject* convert(const std::map<K,V>& map) {
        boost::python::dict result;
        for(typename std::map<K,V>::const_iterator iter=map.begin(); iter!=map.end(); ++iter) {
            result[boost::python::object(iter->first)]=boost::python::object(iter->second);
        }
        return boost::python::incref(boost::python::dict(result).ptr());
    }
    static const PyTypeObject* get_pytype() { return &PyDict_Type; }
};


template<class T1, class T2, class T3>
struct triple_to_python {
    static PyObject* convert(const tuple<T1,T2,T3>& tup) {
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
struct quadruple_to_python {
    static PyObject* convert(const tuple<T1,T2,T3,T4>& tup) {
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
struct array_to_python_list {
    static PyObject* convert(const array<T>& ary) {
        boost::python::list result;
        for(typename array<T>::const_iterator iter=ary.begin(); iter!=ary.end(); ++iter) {
            result.append(boost::python::object(*iter));
        }
        return boost::python::incref(boost::python::list(result).ptr());
    }
    static const PyTypeObject* get_pytype() { return &PyList_Type; }
};

template<class T>
struct list_to_python_list {
    static PyObject* convert(const List<T>& lst) {
        boost::python::list result;
        for(typename List<T>::const_iterator iter=lst.begin(); iter!=lst.end(); ++iter) {
            result.append(boost::python::object(*iter));
        }
        return boost::python::incref(boost::python::list(result).ptr());
    }
    static const PyTypeObject* get_pytype() { return &PyList_Type; }
};

template<class T>
struct set_to_python_list {
    static PyObject* convert(const Set<T>& set) {
        boost::python::list result;
        for(typename Set<T>::const_iterator iter=set.begin(); iter!=set.end(); ++iter) {
            result.append(boost::python::object(*iter));
        }
        return boost::python::incref(boost::python::list(result).ptr());
    }
    static const PyTypeObject* get_pytype() { return &PyList_Type; }
};

template<class T> class Vector;

template<class T>
struct vector_to_python_list {
    static PyObject* convert(const Vector<T>& vec) {
        boost::python::list result;
        for(uint i=0; i!=vec.size(); ++i) {
            result.append(boost::python::object(vec[i]));
        }
        return boost::python::incref(boost::python::list(result).ptr());
    }
    static const PyTypeObject* get_pytype() { return &PyList_Type; }
};





// Read a array variable of type X from a Python object
template<class X>
void
read_list_array(array<X>& ary, const boost::python::object& obj)
{
    // See "Extracting C++ objects" in the Boost Python tutorial
    boost::python::list elements=boost::python::extract<boost::python::list>(obj);
    int n=boost::python::len(elements);
    ary.resize(n);
    for(int i=0; i!=n; ++i) {
        read(ary[i], elements[i]);
    }
}

// Read a array variable of type X from a Python object
template<class X>
void
read_tuple_array(array<X>& ary, const boost::python::object& obj)
{
    // See "Extracting C++ objects" in the Boost Python tutorial
    boost::python::tuple elements=boost::python::extract<boost::python::tuple>(obj);
    int n=boost::python::len(elements);
    ary.resize(n);
    for(int i=0; i!=n; ++i) {
        read(ary[i], elements[i]);
    }
}


template<class T>
T*
make(const boost::python::object& obj)
{
    T* t=new T;
    read(*t,obj);
    return t;
}

template<class T>
T*
make2(const boost::python::object& obj1,const boost::python::object& obj2)
{
    T* t=new T;
    read(*t,obj1,obj2);
    return t;
}

template<class T>
T*
make3(const boost::python::object& obj1,const boost::python::object& obj2,const boost::python::object& obj3)
{
    T* t=new T;
    read(*t,obj1,obj2,obj3);
    return t;
}

template<class C, class I0, class I1, class X> inline
X get_slice(const C& c, const I0& i0, const I1& i1) { return project(c,range(i0,i1)); }

template<class C, class I, class X> inline
X get_item(const C& c, const I& i) { return c[i]; }

template<class C, class I, class X> inline
void set_item(C& c, const I& i, const X& x) { c[i]=x; }

template<class C, class I, class J, class X> inline
X matrix_get_item(const C& c, const I& i, const J& j) { return c[i][j]; }

template<class C, class I, class J, class X> inline
void matrix_set_item(C& c, const I& i, const J& j, const X& x) { c[i][j]=x; }



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
R __prod__(const A1& a1, const A2& a2) { return static_cast<R>(prod(a1,a2)); }

template<class R, class A1, class A2>
R __radd__(const A1& a1, const A2& a2) { return static_cast<R>(a2+a1); }

template<class R, class A1, class A2>
R __rsub__(const A1& a1, const A2& a2) { return static_cast<R>(a2-a1); }

template<class R, class A1, class A2>
R __rmul__(const A1& a1, const A2& a2) { return static_cast<R>(a2*a1); }

template<class R, class A1, class A2>
R __rdiv__(const A1& a1, const A2& a2) { return static_cast<R>(a2/a1); }

/*
template<class R, class A, class T=A>
R __neg__(const A& a) { return static_cast<R>(-static_cast<T>(a)); }

template<class R, class A1, class A2, class T1=A1, class T2=A2>
R __add__(const A1& a1, const A2& a2) { return static_cast<R>(static_cast<T1>(a1)+static_cast<T2>(a2)); }

template<class R, class A1, class A2, class T1=A1, class T2=A2>
R __sub__(const A1& a1, const A2& a2) { return static_cast<R>(static_cast<T1>(a1)-static_cast<T2>(a2)); }

template<class R, class A1, class A2, class T1=A1, class T2=A2>
R __mul__(const A1& a1, const A2& a2) { return static_cast<R>(static_cast<T1>(a1)*static_cast<T2>(a2)); }

template<class R, class A1, class A2, class T1=A1, class T2=A2>
R __div__(const A1& a1, const A2& a2) { return static_cast<R>(static_cast<T1>(a1)/static_cast<T2>(a2)); }

template<class R, class A1, class A2, class T1=A1, class T2=A2>
R __prod__(const A1& a1, const A2& a2) { return static_cast<R>(prod(static_cast<T1>(a1),static_cast<T2>(a2))); }
*/

template<class T> bool check(const boost::python::extract<T>& e) { e(); return e.check(); }
template<> bool check(const boost::python::extract<boost::python::list>& e);
template<> bool check(const boost::python::extract<boost::python::dict>& e);
template<> bool check(const boost::python::extract<boost::python::tuple>& e);


}

#endif /* ARIADNE_PYTHON_UTILITIES_H */
