/***************************************************************************
 *            utilities.hpp
 *
 *  Copyright  2005-20  Alberto Casagrande, Pieter Collins
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

#ifndef ARIADNE_PYTHON_ARIADNE_UTILITIES_HPP
#define ARIADNE_PYTHON_ARIADNE_UTILITIES_HPP

#include "pybind11.hpp"
#include "interval-utilities.hpp"

#include "utility/array.hpp"
#include "utility/tuple.hpp"
#include "utility/container.hpp"
#include "algebra/declarations.hpp"
#include "function/declarations.hpp"
#include "geometry/declarations.hpp"
#include "utility/metaprogramming.hpp"



namespace Ariadne {

using namespace PyBind11;


inline uint pyindex(int i, uint n) { return (i>=0) ? static_cast<uint>(i) : (n-static_cast<uint>(-i)); }

template<class C, class I, class X=decltype(declval<const C>()[declval<I>()])> inline
auto __getitem__(const C& c, const I& i) -> X {
    if constexpr (std::same_as<I,int>) { return c[pyindex(i,c.size())]; } else { return c[i]; } }

template<class C, class I0, class I1, class S> inline
S __getslice__(const C& c, const I0& i0, const I1& i1) {
        if constexpr (std::same_as<I0,int>) { auto n=c.size(); return project(c,range(pyindex(i0,n,pyindex(i1,n)))); } else { return project(c,range(i0,i1)); } }

template<class C, class I, class X> inline
Void __setitem__(C& c, const I& i, const X& x) {
        if constexpr (std::same_as<I,int>) { c[pyindex(i,c.size())]=x; } else { c[i]=x; } }


template<class C, class I, class J, class X> inline
X __getitem2__(const C& c, const I& i, const J& j) { return c[i][j]; }

template<class C, class I, class J, class X> inline
Void __setitem2__(C& c, const I& i, const J& j, const X& x) { c[i][j]=x; }



template<class A1, class A2>
A1& __iadd__(A1& a1, const A2& a2) { return a1+=a2; }

template<class A1, class A2>
A1& __isub__(A1& a1, const A2& a2) { return a1-=a2; }

template<class A1, class A2>
A1& __imul__(A1& a1, const A2& a2) { return a1*=a2; }

template<class A1, class A2>
A1& __idiv__(A1& a1, const A2& a2) { return a1/=a2; }

template<class F, class... AS> auto __call__(F const& f, AS... as) -> decltype(f(as...)) { return f(as...); }
template<class... AS> auto _evaluate_(AS... as) -> decltype(evaluate(as...)) { return evaluate(as...); }
template<class... AS> auto _partial_evaluate_(AS... as) -> decltype(partial_evaluate(as...)) { return partial_evaluate(as...); }
template<class... AS> auto _unchecked_evaluate_(AS... as) -> decltype(unchecked_evaluate(as...)) { return unchecked_evaluate(as...); }
template<class... AS> auto _compose_(AS... as) -> decltype(compose(as...)) { return compose(as...); }
template<class... AS> auto _unchecked_compose_(AS... as) -> decltype(unchecked_compose(as...)) { return unchecked_compose(as...); }

template<class... AS> auto _inverse_(AS... as) -> decltype(inverse(as...)) { return inverse(as...); }
template<class... AS> auto _solve_(AS... as) -> decltype(solve(as...)) { return solve(as...); }

template<class... AS> auto _differential_(AS... as) -> decltype(differential(as...)) { return differential(as...); }

template<class... AS> auto _norm_(AS... as) -> decltype(norm(as...)) { return norm(as...); }
template<class... AS> auto _sup_norm_(AS... as) -> decltype(sup_norm(as...)) { return sup_norm(as...); }
template<class... AS> auto _two_norm_(AS... as) -> decltype(two_norm(as...)) { return two_norm(as...); }
template<class... AS> auto _dot_(AS... as) -> decltype(dot(as...)) { return dot(as...); }
template<class... AS> auto _join_(AS... as) -> decltype(join(as...)) { return join(as...); }
template<class... AS> auto _cojoin_(AS... as) -> decltype(cojoin(as...)) { return cojoin(as...); }
template<class... AS> auto _combine_(AS... as) -> decltype(combine(as...)) { return combine(as...); }
template<class... AS> auto _transpose_(AS... as) -> decltype(transpose(as...)) { return transpose(as...); }

template<class... AS> auto _midpoint_(AS... as) -> decltype(midpoint(as...)) { return midpoint(as...); }
template<class... AS> auto _embed_(AS... as) -> decltype(embed(as...)) { return embed(as...); }
template<class... AS> auto _extension_(AS... as) -> decltype(extension(as...)) { return extension(as...); }
template<class... AS> auto _restriction_(AS... as) -> decltype(restriction(as...)) { return restriction(as...); }
template<class... AS> auto _split_(AS... as) -> decltype(split(as...)) { return split(as...); }
template<class... AS> auto _derivative_(AS... as) -> decltype(derivative(as...)) { return derivative(as...); }
template<class... AS> auto _antiderivative_(AS... as) -> decltype(antiderivative(as...)) { return antiderivative(as...); }

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


} // namespace Ariadnelist


template<class T>
void export_array(pybind11::module& module, const char* name)
{
    using namespace Ariadne;

    pybind11::class_<Array<T>> array_class(module,name);
    array_class.def(pybind11::init<Array<T>>());
    if constexpr (DefaultConstructible<T>) {
        array_class.def(pybind11::init<uint>());
    }
    array_class.def(pybind11::init<uint,T>());
    array_class.def("__len__", &Array<T>::size);
    array_class.def("__getitem__", &__getitem__<Array<T>,int,T>);
    array_class.def("__setitem__", &__setitem__<Array<T>,int,T>);
    array_class.def("__str__", &__cstr__<Array<T>>);
}


namespace pybind11::detail {

// The third template argument is 'true' if the array is resizable
template <class T> struct type_caster<Ariadne::Array<T>>
    : array_caster<Ariadne::Array<T>, T, Ariadne::DefaultConstructible<T>> { };
template <class T> struct type_caster<Ariadne::List<T>>
    : list_caster<Ariadne::List<T>, T> { };
template <class T> struct type_caster<Ariadne::Set<T>>
    : set_caster<Ariadne::Set<T>, T> { };
template <class K, class V> struct type_caster<Ariadne::Map<K,V>>
    : map_caster<Ariadne::Map<K,V>, K,V> { };
} // namespace pybind11::detail



namespace Ariadne {

template<class F, class T> void define_conversion(pybind11::class_<T>& pyclass) {
    if constexpr (Constructible<T,F> and not Same<T,F>) {
        pyclass.def(pybind11::init<F>());
        if constexpr (Convertible<F,T>) {
            pybind11::implicitly_convertible<F,T>();
        }
    }
}

template<class X, class Y>
pybind11::class_<X>& define_inplace_arithmetic(pybind11::module& module, pybind11::class_<X>& pyclass, Tag<Y> = Tag<Y>()) {
    module.def("__iadd__", &__iadd__<X,Y>);
    module.def("__isub__", &__isub__<X,Y>);
    module.def("__imul__", &__imul__<X,Y>);
    module.def("__idiv__", &__idiv__<X,Y>);
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
    if constexpr(HasEquality<X,X>) {
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
pybind11::class_<V>& define_vector_arithmetic(pybind11::module&, pybind11::class_<V>& pyclass, Tag<X> = Tag<X>()) {
    pyclass.def("__pos__", &__pos__<V>, pybind11::is_operator());
    pyclass.def("__neg__", &__neg__<V>, pybind11::is_operator());
    pyclass.def("__add__", &__add__<V,V>, pybind11::is_operator());
    pyclass.def("__radd__", &__radd__<V,V>, pybind11::is_operator());
    pyclass.def("__sub__", &__sub__<V,V>, pybind11::is_operator());
    pyclass.def("__rsub__", &__rsub__<V,V>, pybind11::is_operator());
    pyclass.def("__rmul__", &__rmul__<V,X>, pybind11::is_operator());
    pyclass.def("__mul__", &__mul__<V,X>, pybind11::is_operator());
    if constexpr(CanDivide<X,X>) {
        pyclass.def(__py_div__,__div__<V,X>, pybind11::is_operator());
    }
//    module.def("dot",  &_dot_<Vector<X>,Vector<X>>);
    return pyclass;
}

template<class VX, class VY>
pybind11::class_<VX>& define_mixed_vector_arithmetic(pybind11::module& module, pybind11::class_<VX>& pyclass, Tag<VY> = Tag<VY>()) {
    using X=typename VX::ScalarType;
    using Y=typename VY::ScalarType;
    pyclass.def("__add__", &__add__<VX,VY>, pybind11::is_operator());
    pyclass.def("__radd__", &__radd__<VX,VY>, pybind11::is_operator());
    pyclass.def("__sub__", &__sub__<VX,VY>, pybind11::is_operator());
    pyclass.def("__rsub__", &__rsub__<VX,VY>, pybind11::is_operator());
    pyclass.def("__rmul__", &__rmul__<VX,Y>, pybind11::is_operator());
    pyclass.def("__mul__", &__mul__<VX,Y>, pybind11::is_operator());
    if constexpr(CanDivide<X,Y>) {
        pyclass.def(__py_div__,__div__<VX,Y>, pybind11::is_operator());
    }
    module.def("dot",  &_dot_<Vector<X>,Vector<Y>>);
    module.def("dot",  &_dot_<Vector<Y>,Vector<X>>);
    return pyclass;
}

template<class V, class X=typename V::ScalarType>
pybind11::class_<V>& define_inplace_vector_arithmetic(pybind11::module& module, pybind11::class_<V>& pyclass, Tag<X> = Tag<X>()) {
    module.def("__iadd__", &__iadd__<V,V>, pybind11::is_operator());
    module.def("__isub__", &__isub__<V,V>, pybind11::is_operator());
    module.def("__imul__", &__imul__<V,X>, pybind11::is_operator());
    module.def("__idiv__", &__idiv__<V,X>, pybind11::is_operator());
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
pybind11::class_<VA>& define_vector_algebra_arithmetic(pybind11::module&, pybind11::class_<VA>& pyclass) {
    typedef Vector<X> VX;

    pyclass.def("__pos__", &__pos__<VA>, pybind11::is_operator());
    pyclass.def("__neg__", &__neg__<VA>, pybind11::is_operator());
    pyclass.def("__add__", &__add__<VA,VA>, pybind11::is_operator());
    pyclass.def("__sub__", &__sub__<VA,VA>, pybind11::is_operator());

    pyclass.def("__add__", &__add__<VA,VX>, pybind11::is_operator());
    pyclass.def("__sub__", &__add__<VA,VX>, pybind11::is_operator());
    pyclass.def("__radd__", &__radd__<VA,VX>, pybind11::is_operator());
    pyclass.def("__rsub__", &__radd__<VA,VX>, pybind11::is_operator());

    pyclass.def("__rmul__", &__rmul__<VA,A>, pybind11::is_operator());
    pyclass.def("__mul__", &__mul__<VA,A>, pybind11::is_operator());
    if constexpr(CanDivide<VA,A>) {
        pyclass.def(__py_div__, &__div__<VA,A>, pybind11::is_operator());
    }

    pyclass.def("__rmul__", &__rmul__<VA,X>, pybind11::is_operator());
    pyclass.def("__mul__", &__mul__<VA,X>, pybind11::is_operator());
    if constexpr(CanDivide<VA,X>) {
        pyclass.def(__py_div__, &__div__<VA,X>, pybind11::is_operator());
    }
    return pyclass;
}


} // namespace Ariadne


#include "algebra/vector.hpp"

template<class X>
pybind11::class_<Ariadne::Vector<X>> export_vector(pybind11::module& module, std::string name) {
    using namespace Ariadne;

    pybind11::class_<Vector<X>> vector_class(module, name.c_str());
    vector_class.def(pybind11::init<Vector<X>>());
//    vector_class.def(pybind11::init<Array<X>>());
    if constexpr (DefaultConstructible<X>) {
        vector_class.def(pybind11::init<Nat>());
    }
    vector_class.def(pybind11::init<Nat,X>());
    vector_class.def("size", &Vector<X>::size);
    vector_class.def("__len__", &Vector<X>::size);
    vector_class.def("__setitem__", &__setitem__<Vector<X>,Nat,X>);
    vector_class.def("__getitem__", &__getitem__<Vector<X>,Nat,X>);
    //vector_class.def("__getslice__", &__getslice__<Vector<X>,int,int,Vector<X>>);
    if constexpr(HasEquality<X,X>) {
        vector_class.def("__eq__", &__eq__<Vector<X>,Vector<X> , Return<EqualityType<X,X>> >);
        vector_class.def("__ne__", &__ne__<Vector<X>,Vector<X> , Return<InequalityType<X,X>> >);
    }
    vector_class.def("__pos__", &__pos__<Vector<X> , Return<Vector<X>> >, pybind11::is_operator());
    vector_class.def("__neg__", &__neg__<Vector<X> , Return<Vector<X>> >, pybind11::is_operator());
    vector_class.def("__add__",__add__<Vector<X>,Vector<X> , Return<Vector<SumType<X,X>>> >, pybind11::is_operator());
    vector_class.def("__sub__",__sub__<Vector<X>,Vector<X> , Return<Vector<DifferenceType<X,X>>> >, pybind11::is_operator());
    vector_class.def("__rmul__",__rmul__<Vector<X>,X , Return<Vector<ProductType<X,X>>> >, pybind11::is_operator());
    vector_class.def("__mul__",__mul__<Vector<X>,X , Return<Vector<ProductType<X,X>>> >, pybind11::is_operator());
    if constexpr(CanDivide<X,X>) {
        vector_class.def(__py_div__,__div__<Vector<X>,X , Return<Vector<QuotientType<X,X>>> >, pybind11::is_operator());
    }
    vector_class.def("__str__",&__cstr__<Vector<X>>);
    //vector_class.def("__repr__",&__repr__<Vector<X>>);
    if constexpr (DefaultConstructible<X>) {
        vector_class.def_static("unit",(Vector<X>(*)(SizeType,SizeType))&Vector<X>::unit);
        vector_class.def_static("basis",(Array<Vector<X>>(*)(SizeType))&Vector<X>::basis);
    }

    module.def("dot", &_dot_<Vector<X>,Vector<X>>);

    module.def("join", &_join_<Vector<X>,Vector<X>>);
    module.def("join", &_join_<Vector<X>,X>);
    module.def("join", &_join_<X,Vector<X>>);
    module.def("join", [](X const& x1, X const& x2){return Vector<X>({x1,x2});});

    return vector_class;
}


#endif /* ARIADNE_PYTHON_ARIADNE_UTILITIES_HPP */
