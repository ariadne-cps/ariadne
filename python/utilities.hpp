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
#include "utility/declarations.hpp"




namespace Ariadne {

static const char* __py_div__ = (PY_MAJOR_VERSION>=3) ? "__truediv__" : "__div__";
static const char* __py_rdiv__ = (PY_MAJOR_VERSION>=3) ? "__rtruediv__" : "__rdiv__";

class String;
template<class X> String class_name();

template<class X> inline std::string python_name(const std::string& name) {
    return class_name<X>()+name; }

inline uint pyindex(int i, uint n) { return (i>=0) ? static_cast<uint>(i) : (n-static_cast<uint>(-i)); }

template<class C, class I, class X=decltype(declval<const C>()[declval<I>()])> inline
auto __getitem__(const C& c, const I& i) -> X {
    if constexpr (std::is_same<I,int>::value) { return c[pyindex(i,c.size())]; } else { return c[i]; } }

template<class C, class I0, class I1, class S> inline
S __getslice__(const C& c, const I0& i0, const I1& i1) {
        if constexpr (std::is_same<I0,int>::value) { auto n=c.size(); return project(c,range(pyindex(i0,n,pyindex(i1,n)))); } else { return project(c,range(i0,i1)); } }

template<class C, class I, class X> inline
Void __setitem__(C& c, const I& i, const X& x) {
        if constexpr (std::is_same<I,int>::value) { c[pyindex(i,c.size())]=x; } else { c[i]=x; } }


template<class C, class I, class J, class X> inline
X __getitem2__(const C& c, const I& i, const J& j) { return c[i][j]; }

template<class C, class I, class J, class X> inline
Void __setitem2__(C& c, const I& i, const J& j, const X& x) { c[i][j]=x; }



template<class A, class RET=Return<decltype(+declval<A>())>>
ReturnType<RET> __pos__(const A& a) { return static_cast<ReturnType<RET>>(+a); }

template<class A, class RET=Return<decltype(-declval<A>())>>
ReturnType<RET> __neg__(const A& a) { return static_cast<ReturnType<RET>>(-a); }

template<class A1, class A2, class RET=Return<decltype(declval<A1>()+declval<A2>())>>
ReturnType<RET> __add__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(a1+a2); }

template<class A1, class A2, class RET=Return<decltype(declval<A1>()-declval<A2>())>>
ReturnType<RET> __sub__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(a1-a2); }

template<class A1, class A2, class RET=Return<decltype(declval<A1>()*declval<A2>())>>
ReturnType<RET> __mul__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(a1*a2); }

template<class A1, class A2, class RET=Return<decltype(declval<A1>()/declval<A2>())>>
ReturnType<RET> __div__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(a1/a2); }

template<class A1, class A2, class RET=Return<decltype(declval<A2>()+declval<A1>())>>
ReturnType<RET> __radd__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(a2+a1); }

template<class A1, class A2, class RET=Return<decltype(declval<A2>()-declval<A1>())>>
ReturnType<RET> __rsub__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(a2-a1); }

template<class A1, class A2, class RET=Return<decltype(declval<A2>()*declval<A1>())>>
ReturnType<RET> __rmul__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(a2*a1); }

template<class A1, class A2, class RET=Return<decltype(declval<A2>()/declval<A1>())>>
ReturnType<RET> __rdiv__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(a2/a1); }

template<class A1, class A2, class RET=Return<decltype(pow(declval<A1>(),declval<A2>()))>>
ReturnType<RET> __pow__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(pow(a1,a2)); }

template<class A1, class A2, class RET=Return<decltype(declval<A1>()&&declval<A2>())>>
ReturnType<RET> __and__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(a1 && a2); }

template<class A1, class A2, class RET=Return<decltype(declval<A1>()||declval<A2>())>>
ReturnType<RET> __or__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(a1 || a2); }

template<class A1, class A2, class RET=Return<decltype(declval<A1>()|declval<A2>())>>
ReturnType<RET> __bitor__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(a1 | a2); }

template<class A, class RET=Return<decltype(!declval<A>())>>
ReturnType<RET> __not__(const A& a) { return static_cast<ReturnType<RET>>(!a); }

template<class A1, class A2, class RET=Return<decltype(declval<A1>()==declval<A2>())>>
ReturnType<RET> __eq__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(a1==a2); }

template<class A1, class A2, class RET=Return<decltype(declval<A1>()!=declval<A2>())>>
ReturnType<RET> __ne__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(a1!=a2); }

template<class A1, class A2, class RET=Return<decltype(declval<A1>()>declval<A2>())>>
ReturnType<RET> __gt__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(a1>a2); }

template<class A1, class A2, class RET=Return<decltype(declval<A1>()<declval<A2>())>>
ReturnType<RET> __lt__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(a1<a2); }

template<class A1, class A2, class RET=Return<decltype(declval<A1>()>=declval<A2>())>>
ReturnType<RET> __ge__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(a1>=a2); }

template<class A1, class A2, class RET=Return<decltype(declval<A1>()<=declval<A2>())>>
ReturnType<RET> __le__(const A1& a1, const A2& a2) { return static_cast<ReturnType<RET>>(a1<=a2); }


template<class A1, class A2>
A1& __iadd__(A1& a1, const A2& a2) { return a1+=a2; }

template<class A1, class A2>
A1& __isub__(A1& a1, const A2& a2) { return a1-=a2; }

template<class A1, class A2>
A1& __imul__(A1& a1, const A2& a2) { return a1*=a2; }

template<class A1, class A2>
A1& __idiv__(A1& a1, const A2& a2) { return a1/=a2; }

template<class... AS> auto _add_(AS const& ... as) -> decltype(add(as...)) { return add(as...); }
template<class... AS> auto _sub_(AS const& ... as) -> decltype(sub(as...)) { return sub(as...); }
template<class... AS> auto _mul_(AS const& ... as) -> decltype(mul(as...)) { return mul(as...); }
template<class... AS> auto _div_(AS const& ... as) -> decltype(div(as...)) { return div(as...); }
template<class... AS> auto _fma_(AS const& ... as) -> decltype(fma(as...)) { return fma(as...); }

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

template<class A> auto _log2_(A const& a) -> decltype(log2(a)) { return log2(a); }

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
    array_class.def("__getitem__", &__getitem__<Array<T>,int,T>);
    array_class.def("__setitem__", &__setitem__<Array<T>,int,T>);
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



namespace Ariadne {

template<class... TS> struct Tag { };

template<class A> pybind11::class_<A>& define_logical(pybind11::module& module, pybind11::class_<A>& pyclass) {
    pyclass.def("__and__", &__and__<A,A>);
    pyclass.def("__or__", &__or__<A,A>);
    pyclass.def("__not__", &__not__<A>);
    return pyclass;
}

template<class A> pybind11::class_<A>& define_lattice(pybind11::module& module, pybind11::class_<A>& pyclass) {
    module.def("abs", &_abs_<A>);
    module.def("max", &_max_<A,A>);
    module.def("min", &_min_<A,A>);
    return pyclass;
}

template<class A> pybind11::class_<A>& define_comparisons(pybind11::module& module, pybind11::class_<A>& pyclass) {
    pyclass.def("__eq__", &__eq__<A,A>);
    pyclass.def("__ne__", &__ne__<A,A>);
    pyclass.def("__le__", &__le__<A,A>);
    pyclass.def("__ge__", &__ge__<A,A>);
    pyclass.def("__lt__", &__lt__<A,A>);
    pyclass.def("__gt__", &__lt__<A,A>);
    return pyclass;
}


template<class X> pybind11::class_<X>& define_arithmetic(pybind11::module& module, pybind11::class_<X>& pyclass) {
    pyclass.def("__pos__", &__pos__<X>);
    pyclass.def("__neg__", &__neg__<X>);
    pyclass.def("__add__", &__add__<X,X>);
    pyclass.def("__sub__", &__sub__<X,X>);
    pyclass.def("__mul__", &__mul__<X,X>);
    if constexpr(CanDivide<X,X>::value) {
        pyclass.def(__py_div__, &__div__<X,X>);
    }
    return pyclass;
}

template<class X, class Y> pybind11::class_<X>& define_mixed_arithmetic(pybind11::module& module, pybind11::class_<X>& pyclass, Tag<Y> = Tag<Y>()) {
    pyclass.def("__add__", &__add__<X,Y>);
    pyclass.def("__radd__", &__radd__<X,Y>);
    pyclass.def("__sub__", &__sub__<X,Y>);
    pyclass.def("__rsub__", &__rsub__<X,Y>);
    pyclass.def("__mul__", &__mul__<X,Y>);
    pyclass.def("__rmul__", &__rmul__<X,Y>);
    if constexpr(CanDivide<X,Y>::value) {
        pyclass.def(__py_div__, &__div__<X,Y>);
    }
    if constexpr(CanDivide<Y,X>::value) {
        pyclass.def("__py_rdiv__", &__rdiv__<X,Y>);
    }
    return pyclass;
}

template<class X, class Y>
pybind11::class_<X>& define_inplace_arithmetic(pybind11::module& module, pybind11::class_<X>& pyclass, Tag<Y> = Tag<Y>()) {
    module.def("__iadd__", &__iadd__<X,Y>);
    module.def("__isub__", &__isub__<X,Y>);
    module.def("__imul__", &__imul__<X,Y>);
    module.def("__idiv__", &__idiv__<X,Y>);
    return pyclass;
}


template<class X> pybind11::class_<X>& define_transcendental(pybind11::module& module, pybind11::class_<X>& pyclass) {
    module.def("sqr", &_sqr_<X>);
    module.def("rec", &_rec_<X>);
    module.def("pow", &_pow_<X,Int>);
    module.def("sqrt", &_sqrt_<X>);
    module.def("exp", &_exp_<X>);
    module.def("log", &_log_<X>);
    module.def("sin", &_sin_<X>);
    module.def("cos", &_cos_<X>);
    module.def("tan", &_tan_<X>);
    module.def("atan", &_atan_<X>);
    return pyclass;
}

template<class X> pybind11::class_<X>& define_elementary(pybind11::module& module, pybind11::class_<X>& pyclass) {
    define_arithmetic(module,pyclass);
    define_transcendental(module,pyclass);
    return pyclass;
}


template<class X> pybind11::class_<X>& define_monotonic(pybind11::module& module, pybind11::class_<X>& pyclass) {
    using NX = decltype(-declval<X>());
    pyclass.def("__pos__", &__pos__<X>);
    pyclass.def("__neg__", &__neg__<X>);
    pyclass.def("__add__", &__add__<X,X>);
    pyclass.def("__sub__", &__sub__<X,NX>);
    module.def("sqrt", &_sqrt_<X>);
    module.def("exp", &_exp_<X>);
    module.def("log", &_log_<X>);
    module.def("atan", &_atan_<X>);
    return pyclass;
}

template<class X, class Y> pybind11::class_<X>& define_mixed_monotonic(pybind11::module& module, pybind11::class_<X>& pyclass, Tag<Y> = Tag<Y>()) {
    using NY = decltype(-declval<Y>());
    pyclass.def("__add__", &__add__<X,Y>);
    pyclass.def("__radd__", &__radd__<X,Y>);
    pyclass.def("__sub__", &__sub__<X,NY>);
    pyclass.def("__rsub__", &__rsub__<X,NY>);
    return pyclass;
}


template<class A, class X=typename A::NumericType> pybind11::class_<A>& define_algebra(pybind11::module& module, pybind11::class_<A>& pyclass, Tag<X> = Tag<X>()) {
    define_arithmetic(module,pyclass);
    define_mixed_arithmetic(module,pyclass,Tag<X>());
    return pyclass;
}


template<class A, class X=typename A::NumericType> pybind11::class_<A>& define_elementary_algebra(pybind11::module& module, pybind11::class_<A>& pyclass, Tag<X> = Tag<X>()) {
    define_algebra<A,X>(module,pyclass);
    define_transcendental(module,pyclass);
    return pyclass;
}

template<class A, class X=typename A::NumericType>
pybind11::class_<A>& define_inplace_algebra(pybind11::module& module, pybind11::class_<A>& pyclass, Tag<X> = Tag<X>()) {
    module.def("__iadd__", &__iadd__<A,A>);
    module.def("__iadd__", &__isub__<A,A>);
    define_inplace_arithmetic(module,pyclass,Tag<X>());
    return pyclass;
}


template<class V, class X=typename V::ScalarType>
pybind11::class_<V>& define_vector_operations(pybind11::module& module, pybind11::class_<V>& pyclass)
{
    pyclass.def("size", &V::size);
    pyclass.def("__len__", &V::size);
    pyclass.def("__setitem__", &__setitem__<V,Int,X>);
    pyclass.def("__getitem__", &__getitem__<V,Int>);
    if constexpr(HasEquality<X,X>::value) {
        pyclass.def("__eq__", &__eq__<V,V , Return<EqualityType<X,X>> >);
        pyclass.def("__ne__", &__ne__<V,V , Return<InequalityType<X,X>> >); }
    pyclass.def("__str__",&__cstr__<V>);
    pyclass.def("__repr__",&__repr__<V>);

    module.def("join", &_join_<V,V>);
    module.def("join", &_join_<V,X>);
    module.def("join", &_join_<X,V>);

//    module.def("join", &__sjoin__<X,X>);

    return pyclass;
}

template<class V, class X=typename V::ScalarType>
pybind11::class_<V>& define_vector_arithmetic(pybind11::module& module, pybind11::class_<V>& pyclass, Tag<X> = Tag<X>()) {
    pyclass.def("__pos__", &__pos__<V>);
    pyclass.def("__neg__", &__neg__<V>);
    pyclass.def("__add__", &__add__<V,V>);
    pyclass.def("__sub__", &__sub__<V,V>);
    pyclass.def("__rmul__", &__rmul__<V,X>);
    pyclass.def("__mul__", &__mul__<V,X>);
    if constexpr(CanDivide<X,X>::value) {
        pyclass.def(__py_div__,__div__<V,X>);
    }
//    module.def("dot",  &_dot_<Vector<X>,Vector<X>>);
    return pyclass;
}

template<class VX, class VY>
pybind11::class_<VX>& define_mixed_vector_arithmetic(pybind11::module& module, pybind11::class_<VX>& pyclass, Tag<VY> = Tag<VY>()) {
    using X=typename VX::ScalarType;
    using Y=typename VY::ScalarType;
    pyclass.def("__add__", &__add__<VX,VY>);
    pyclass.def("__radd__", &__radd__<VX,VY>);
    pyclass.def("__sub__", &__sub__<VX,VY>);
    pyclass.def("__rsub__", &__rsub__<VX,VY>);
    pyclass.def("__rmul__", &__rmul__<VX,Y>);
    pyclass.def("__mul__", &__mul__<VX,Y>);
    if constexpr(CanDivide<X,Y>::value) {
        pyclass.def(__py_div__,__div__<VX,Y>);
    }
    module.def("dot",  &_dot_<Vector<X>,Vector<Y>>);
    module.def("dot",  &_dot_<Vector<Y>,Vector<X>>);
    return pyclass;
}

template<class V, class X=typename V::ScalarType>
pybind11::class_<V>& define_inplace_vector_arithmetic(pybind11::module& module, pybind11::class_<V>& pyclass, Tag<X> = Tag<X>()) {
    module.def("__iadd__", &__iadd__<V,V>);
    module.def("__isub__", &__isub__<V,V>);
    module.def("__imul__", &__imul__<V,X>);
    module.def("__idiv__", &__idiv__<V,X>);
    return pyclass;
}

template<class VX, class VY>
pybind11::class_<VX>& define_inplace_mixed_vector_arithmetic(pybind11::module& module, pybind11::class_<VX>& pyclass, Tag<VY> = Tag<VY>()) {
    using Y=typename VY::ScalarType;
    module.def("__iadd__", &__iadd__<VX,VY>);
    module.def("__isub__", &__isub__<VX,VY>);
    module.def("__imul__", &__imul__<VX,Y>);
    module.def("__idiv__", &__idiv__<VX,Y>);
    return pyclass;
}

template<class V, class X=typename V::ScalarType>
pybind11::class_<V>& define_vector_concept(pybind11::module& module, pybind11::class_<V>& pyclass)
{
    define_vector_operations<V,X>(module,pyclass);
    define_vector_arithmetic<V,X>(module,pyclass);
    return pyclass;
}


template<class VA, class A=typename VA::ScalarType, class X=typename A::NumericType>
pybind11::class_<VA>& define_vector_algebra_arithmetic(pybind11::module& module, pybind11::class_<VA>& pyclass) {
    typedef Vector<X> VX;

    pyclass.def("__pos__", &__pos__<VA>);
    pyclass.def("__neg__", &__neg__<VA>);
    pyclass.def("__add__", &__add__<VA,VA>);
    pyclass.def("__sub__", &__sub__<VA,VA>);

    pyclass.def("__add__", &__add__<VA,VX>);
    pyclass.def("__sub__", &__add__<VA,VX>);
    pyclass.def("__radd__", &__radd__<VA,VX>);
    pyclass.def("__rsub__", &__radd__<VA,VX>);

    pyclass.def("__rmul__", &__rmul__<VA,A>);
    pyclass.def("__mul__", &__mul__<VA,A>);
    if constexpr(CanDivide<VA,A>::value) {
        pyclass.def(__py_div__, &__div__<VA,A>);
    }

    pyclass.def("__rmul__", &__rmul__<VA,X>);
    pyclass.def("__mul__", &__mul__<VA,X>);
    if constexpr(CanDivide<VA,X>::value) {
        pyclass.def(__py_div__, &__div__<VA,X>);
    }
    return pyclass;
}


} // namespace Ariadne


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
        vector_class.def("__eq__", &__eq__<Vector<X>,Vector<X> , Return<EqualityType<X,X>> >);
        vector_class.def("__ne__", &__ne__<Vector<X>,Vector<X> , Return<InequalityType<X,X>> >);
    }
    vector_class.def("__pos__", &__pos__<Vector<X> , Return<Vector<X>> >);
    vector_class.def("__neg__", &__neg__<Vector<X> , Return<Vector<X>> >);
    vector_class.def("__add__",__add__<Vector<X>,Vector<X> , Return<Vector<SumType<X,X>>> >);
    vector_class.def("__sub__",__sub__<Vector<X>,Vector<X> , Return<Vector<DifferenceType<X,X>>> >);
    vector_class.def("__rmul__",__rmul__<Vector<X>,X , Return<Vector<ProductType<X,X>>> >);
    vector_class.def("__mul__",__mul__<Vector<X>,X , Return<Vector<ProductType<X,X>>> >);
    if constexpr(CanDivide<X,X>::value) {
        vector_class.def(__py_div__,__div__<Vector<X>,X , Return<Vector<QuotientType<X,X>>> >);
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
