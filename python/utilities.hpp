/***************************************************************************
 *            utilities.hpp
 *
 *  Copyright  2005-8  Alberto Casagrande, Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file utilities.hpp
 *  Commonly used inline methods for the Python interface.
 */

#ifndef ARIADNE_PYTHON_UTILITIES_HPP
#define ARIADNE_PYTHON_UTILITIES_HPP

#include "pybind11.hpp"

#include "utility/array.hpp"
#include "utility/tuple.hpp"
#include "utility/container.hpp"

namespace Ariadne {

class String;
template<class X> String class_name();

template<class X> inline std::string python_name(const std::string& name) {
    return class_name<X>()+name; }



template<class C, class I0, class I1, class S> inline
S __getslice__(const C& c, const I0& i0, const I1& i1) { return project(c,range(i0,i1)); }

template<class C, class I, class X, EnableIf<IsSame<I,int>> =dummy> inline
X __getitem__(const C& c, const I& i) {
    return c[static_cast<uint>(i)];
}

template<class C, class I, class X, EnableIf<Not<IsSame<I,int>>> =dummy> inline
X __getitem__(const C& c, const I& i) {
    return c[i];
}


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


template<class... AS> decltype(auto) _nul_(AS const& ... as) { return nul(as...); }
template<class... AS> decltype(auto) _pos_(AS const& ... as) { return pos(as...); }
template<class... AS> decltype(auto) _neg_(AS const& ... as) { return neg(as...); }
template<class... AS> decltype(auto) _hlf_(AS const& ... as) { return hlf(as...); }
template<class... AS> decltype(auto) _rec_(AS const& ... as) { return rec(as...); }
template<class... AS> decltype(auto) _sqr_(AS const& ... as) { return sqr(as...); }
template<class... AS> decltype(auto) _pow_(AS const& ... as) { return pow(as...); }
template<class... AS> decltype(auto) _sqrt_(AS const& ... as) { return sqrt(as...); }
template<class... AS> decltype(auto) _exp_(AS const& ... as) { return exp(as...); }
template<class... AS> decltype(auto) _log_(AS const& ... as) { return log(as...); }
template<class... AS> decltype(auto) _sin_(AS const& ... as) { return sin(as...); }
template<class... AS> decltype(auto) _cos_(AS const& ... as) { return cos(as...); }
template<class... AS> decltype(auto) _tan_(AS const& ... as) { return tan(as...); }
template<class... AS> decltype(auto) _atan_(AS const& ... as) { return atan(as...); }
template<class... AS> decltype(auto) _max_(AS const& ... as) { return max(as...); }
template<class... AS> decltype(auto) _min_(AS const& ... as) { return min(as...); }
template<class... AS> decltype(auto) _abs_(AS const& ... as) { return abs(as...); }
template<class... AS> decltype(auto) _sgn_(AS const& ... as) { return sgn(as...); }

template<class... AS> inline decltype(auto) _evaluate_(AS... as) { return evaluate(as...); }
template<class... AS> inline decltype(auto) _partial_evaluate_(AS... as) { return partial_evaluate(as...); }
template<class... AS> inline decltype(auto) _unchecked_evaluate_(AS... as) { return unchecked_evaluate(as...); }
template<class... AS> inline decltype(auto) _compose_(AS... as) { return compose(as...); }
template<class... AS> inline decltype(auto) _unchecked_compose_(AS... as) { return unchecked_compose(as...); }

template<class... AS> inline decltype(auto) _dot_(AS... as) { return dot(as...); }
template<class... AS> inline decltype(auto) _join_(AS... as) { return join(as...); }
template<class... AS> inline decltype(auto) _combine_(AS... as) { return combine(as...); }

template<class... AS> inline decltype(auto) _midpoint_(AS... as) { return midpoint(as...); }
template<class... AS> inline decltype(auto) _embed_(AS... as) { return embed(as...); }
template<class... AS> inline decltype(auto) _extension_(AS... as) { return extension(as...); }
template<class... AS> inline decltype(auto) _restriction_(AS... as) { return restriction(as...); }
template<class... AS> inline decltype(auto) _split_(AS... as) { return split(as...); }
template<class... AS> inline decltype(auto) _derivative_(AS... as) { return derivative(as...); }
template<class... AS> inline decltype(auto) _antiderivative_(AS... as) { return antiderivative(as...); }

template<class... AS> inline decltype(auto) _refinement_(AS... as) { return refinement(as...); }
template<class... AS> inline decltype(auto) _refines_(AS... as) { return refines(as...); }
template<class... AS> inline decltype(auto) _inconsistent_(AS... as) { return inconsistent(as...); }

template<class... AS> decltype(auto) _widen_(AS const& ... as) { return widen(as...); }

template<class... AS> decltype(auto) _contains_(AS const& ... as) { return contains(as...); }
template<class... AS> decltype(auto) _intersection_(AS const& ... as) { return intersection(as...); }
template<class... AS> decltype(auto) _disjoint_(AS const& ... as) { return disjoint(as...); }
template<class... AS> decltype(auto) _subset_(AS const& ... as) { return subset(as...); }
template<class... AS> decltype(auto) _product_(AS const& ... as) { return product(as...); }
template<class... AS> decltype(auto) _hull_(AS const& ... as) { return hull(as...); }
template<class... AS> decltype(auto) _separated_(AS const& ... as) { return separated(as...); }
template<class... AS> decltype(auto) _overlap_(AS const& ... as) { return overlap(as...); }
template<class... AS> decltype(auto) _covers_(AS const& ... as) { return covers(as...); }
template<class... AS> decltype(auto) _inside_(AS const& ... as) { return inside(as...); }

template<class... AS> decltype(auto) _image_(AS const& ... as) { return image(as...); }
template<class... AS> decltype(auto) _preimage_(AS const& ... as) { return preimage(as...); }


template<class T> std::string __cstr__(const T& t) {
    std::stringstream ss; ss << t; return ss.str(); }

template<class T> std::string __crepr__(const T& t) {
    std::stringstream ss; ss << representation(t); return ss.str(); }


template<class T> struct PythonRepresentation {
    const T* pointer;
    PythonRepresentation(const T& t) : pointer(&t) { }
    const T& reference() const { return *pointer; }
};

template<class T> PythonRepresentation<T> python_representation(const T& t) {
    return PythonRepresentation<T>(t); }

template<class T> std::string __repr__(const T& t) {
    std::stringstream ss; ss << python_representation(t); return ss.str();}

} // namespace Ariadne


template<class T>
void export_array(pybind11::module& module, const char* name)
{
    using namespace Ariadne;

    pybind11::class_<Array<T>> array_class(module,name);
    array_class.def(pybind11::init<Array<T>>());
    array_class.def(pybind11::init<uint>());
    array_class.def("__len__", &Array<T>::size);
    array_class.def("__getitem__", &__getitem__<Array<T>,uint,T>);
    array_class.def("__setitem__", &__setitem__<Array<T>,uint,T>);
    array_class.def("__str__", &__cstr__<Array<T>>);
}


namespace pybind11::detail {
template <class T> struct type_caster<Ariadne::Array<T>>
    : array_caster<Ariadne::Array<T>, T, true> { };
template <class T> struct type_caster<Ariadne::List<T>>
    : list_caster<Ariadne::List<T>, T> { };
template <class T> struct type_caster<Ariadne::Set<T>>
    : set_caster<Ariadne::Set<T>, T> { };
template <class K, class V> struct type_caster<Ariadne::Map<K,V>>
    : map_caster<Ariadne::Map<K,V>, K,V> { };
} // namespace pybind11::detail



#include "algebra/vector.hpp"

template<class X>
void export_vector(pybind11::module& module, std::string name) {
    using namespace Ariadne;

    pybind11::class_<Vector<X>> vector_class(module, name.c_str());
    vector_class.def(pybind11::init<Vector<X>>());
    vector_class.def(pybind11::init<Array<X>>());
    if constexpr (IsDefaultConstructible<X>::value) {
        vector_class.def(pybind11::init<Nat>());
    }
    vector_class.def(pybind11::init<Nat,X>());
    vector_class.def("size", &Vector<X>::size);
    vector_class.def("__len__", &Vector<X>::size);
    vector_class.def("__setitem__", &__setitem__<Vector<X>,Nat,X>);
    vector_class.def("__getitem__", &__getitem__<Vector<X>,Nat,X>);
    //vector_class.def("__getslice__", &__getslice__<Vector<X>,int,int,Vector<X>>);
    if constexpr(HasEquality<X,X>::value) {
        vector_class.def("__eq__", &__eq__<EqualityType<X,X>,Vector<X>,Vector<X> >);
        vector_class.def("__ne__", &__ne__<InequalityType<X,X>,Vector<X>,Vector<X> >);
    }
    vector_class.def("__pos__", &__pos__<Vector<X>,Vector<X>>);
    vector_class.def("__neg__", &__neg__<Vector<X>,Vector<X>>);
    vector_class.def("__add__",__add__<Vector<SumType<X,X>>, Vector<X>, Vector<X> >);
    vector_class.def("__sub__",__sub__<Vector<DifferenceType<X,X>>, Vector<X>, Vector<X> >);
    vector_class.def("__rmul__",__rmul__<Vector<ProductType<X,X>>, Vector<X>, X >);
    vector_class.def("__mul__",__mul__<Vector<ProductType<X,X>>, Vector<X>, X >);
    if constexpr(CanDivide<X,X>::value) {
        vector_class.def("__div__",__div__<Vector<QuotientType<X,X>>, Vector<X>, X >);
    }
    vector_class.def("__str__",&__cstr__<Vector<X>>);
    //vector_class.def("__repr__",&__repr__<Vector<X>>);
    vector_class.def_static("unit",&Vector<X>::unit);
    vector_class.def_static("basis",&Vector<X>::basis);

    module.def("dot", &_dot_<Vector<X>,Vector<X>>);

    module.def("join", &_join_<Vector<X>,Vector<X>>);
    module.def("join", &_join_<Vector<X>,X>);
    module.def("join", &_join_<X,Vector<X>>);
    module.def("join", [](X const& x1, X const& x2){return Vector<X>({x1,x2});});
};


#endif /* ARIADNE_PYTHON_UTILITIES_HPP */
