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


template<class... AS> auto _nul_(AS const& ... as) -> decltype(nul(as...)) { return nul(as...); }
template<class... AS> auto _pos_(AS const& ... as) -> decltype(pos(as...)) { return pos(as...); }
template<class... AS> auto _neg_(AS const& ... as) -> decltype(neg(as...)) { return neg(as...); }
template<class... AS> auto _hlf_(AS const& ... as) -> decltype(hlf(as...)) { return hlf(as...); }
template<class... AS> auto _rec_(AS const& ... as) -> decltype(rec(as...)) { return rec(as...); }
template<class... AS> auto _sqr_(AS const& ... as) -> decltype(sqr(as...)) { return sqr(as...); }
template<class... AS> auto _pow_(AS const& ... as) -> decltype(pow(as...)) { return pow(as...); }
template<class... AS> auto _sqrt_(AS const& ... as) -> decltype(sqrt(as...)) { return sqrt(as...); }
template<class... AS> auto _exp_(AS const& ... as) -> decltype(exp(as...)) { return exp(as...); }
template<class... AS> auto _log_(AS const& ... as) -> decltype(log(as...)) { return log(as...); }
template<class... AS> auto _sin_(AS const& ... as) -> decltype(sin(as...)) { return sin(as...); }
template<class... AS> auto _cos_(AS const& ... as) -> decltype(cos(as...)) { return cos(as...); }
template<class... AS> auto _tan_(AS const& ... as) -> decltype(tan(as...)) { return tan(as...); }
template<class... AS> auto _atan_(AS const& ... as) -> decltype(atan(as...)) { return atan(as...); }
template<class... AS> auto _max_(AS const& ... as) -> decltype(max(as...)) { return max(as...); }
template<class... AS> auto _min_(AS const& ... as) -> decltype(min(as...)) { return min(as...); }
template<class... AS> auto _abs_(AS const& ... as) -> decltype(abs(as...)) { return abs(as...); }
template<class... AS> auto _sgn_(AS const& ... as) -> decltype(sgn(as...)) { return sgn(as...); }

template<class... AS> auto _evaluate_(AS... as) -> decltype(evaluate(as...)) { return evaluate(as...); }
template<class... AS> auto _partial_evaluate_(AS... as) -> decltype(partial_evaluate(as...)) { return partial_evaluate(as...); }
template<class... AS> auto _unchecked_evaluate_(AS... as) -> decltype(unchecked_evaluate(as...)) { return unchecked_evaluate(as...); }
template<class... AS> auto _compose_(AS... as) -> decltype(compose(as...)) { return compose(as...); }
template<class... AS> auto _unchecked_compose_(AS... as) -> decltype(unchecked_compose(as...)) { return unchecked_compose(as...); }

template<class... AS> auto _dot_(AS... as) -> decltype(dot(as...)) { return dot(as...); }
template<class... AS> auto _join_(AS... as) -> decltype(join(as...)) { return join(as...); }
template<class... AS> auto _combine_(AS... as) -> decltype(combine(as...)) { return combine(as...); }

template<class... AS> auto _midpoint_(AS... as) -> decltype(midpoint(as...)) { return midpoint(as...); }
template<class... AS> auto _embed_(AS... as) -> decltype(embed(as...)) { return embed(as...); }
template<class... AS> auto _extension_(AS... as) -> decltype(extension(as...)) { return extension(as...); }
template<class... AS> auto _restriction_(AS... as) -> decltype(restriction(as...)) { return restriction(as...); }
template<class... AS> auto _split_(AS... as) -> decltype(split(as...)) { return split(as...); }
template<class... AS> auto _derivative_(AS... as) -> decltype(derivative(as...)) { return derivative(as...); }
template<class... AS> auto _antiderivative_(AS... as) -> decltype(antiderivative(as...)) { return antiderivative(as...); }

template<class... AS> auto _refinement_(AS... as) -> decltype(refinement(as...)) { return refinement(as...); }
template<class... AS> auto _refines_(AS... as) -> decltype(refines(as...)) { return refines(as...); }
template<class... AS> auto _inconsistent_(AS... as) -> decltype(inconsistent(as...)) { return inconsistent(as...); }

template<class... AS> auto _widen_(AS const& ... as) -> decltype(widen(as...)) { return widen(as...); }

template<class... AS> auto _contains_(AS const& ... as) -> decltype(contains(as...)) { return contains(as...); }
template<class... AS> auto _intersection_(AS const& ... as) -> decltype(intersection(as...)) { return intersection(as...); }
template<class... AS> auto _disjoint_(AS const& ... as) -> decltype(disjoint(as...)) { return disjoint(as...); }
template<class... AS> auto _subset_(AS const& ... as) -> decltype(subset(as...)) { return subset(as...); }
template<class... AS> auto _product_(AS const& ... as) -> decltype(product(as...)) { return product(as...); }
template<class... AS> auto _hull_(AS const& ... as) -> decltype(hull(as...)) { return hull(as...); }
template<class... AS> auto _separated_(AS const& ... as) -> decltype(separated(as...)) { return separated(as...); }
template<class... AS> auto _overlap_(AS const& ... as) -> decltype(overlap(as...)) { return overlap(as...); }
template<class... AS> auto _covers_(AS const& ... as) -> decltype(covers(as...)) { return covers(as...); }
template<class... AS> auto _inside_(AS const& ... as) -> decltype(inside(as...)) { return inside(as...); }

template<class... AS> auto _image_(AS const& ... as) -> decltype(image(as...)) { return image(as...); }
template<class... AS> auto _preimage_(AS const& ... as) -> decltype(preimage(as...)) { return preimage(as...); }


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
}


#endif /* ARIADNE_PYTHON_UTILITIES_HPP */
