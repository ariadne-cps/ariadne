/***************************************************************************
 *            float_operations.hpp
 *
 *  Copyright 2008-17  Pieter Collins
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

/*! \file float_operations.hpp
 *  \brief Common header for floating-point number operations.
 */

#ifndef ARIADNE_FLOAT_OPERATIONS_HPP
#define ARIADNE_FLOAT_OPERATIONS_HPP

#include "../utility/macros.hpp"

#include "number.decl.hpp"
#include "float.decl.hpp"
#include "floatdp.hpp"
#include "floatmp.hpp"
#include "float-raw.hpp"

namespace Ariadne {

template<class X, class R=X> struct DeclareFloatOperations
    : public DeclareRealOperations<X,R>
    , public DeclareComparisonOperations<X,LessTrait<X>,EqualsTrait<X>>
    , public DeclareMixedFieldOperators<X,GenericTrait<X>,R>
{
    friend OutputStream& operator<<(OutputStream&, X const&);
    friend InputStream& operator>>(InputStream&, X&);
};

template<class X, class NX=OppositeTrait<X>> struct DeclareDirectedFloatOperations
    : DeclareDirectedNumericOperations<X,NX>
    , DeclareMixedDirectedGroupOperators<X,NX,GenericTrait<X>,GenericTrait<NX>>
{
    friend OutputStream& operator<<(OutputStream&, X const&);
    friend InputStream& operator>>(InputStream&, X&);
};

template<class PX, class QPX=OppositeTrait<PX>> struct DeclarePositiveDirectedFloatOperations
    : DeclarePositiveDirectedNumericOperations<PX,QPX>
    , DeclareMixedDirectedSemifieldOperators<PX,QPX,GenericTrait<PX>,GenericTrait<QPX>>
{
    friend OutputStream& operator<<(OutputStream&, PX const&);
    friend InputStream& operator>>(InputStream&, PX&);
};

template<class PX> struct DeclarePositiveFloatOperations
    : DeclarePositiveDirectedFloatOperations<PX,PX>
{
};

template<class X, class R=X> struct DispatchFloatOperations
    : DispatchNumericOperations<X,R>
    , DispatchComparisonOperations<X,LessTrait<X>,EqualsTrait<X>>
    , DefineConcreteGenericOperators<X>
//    , ProvideConcreteGenericFieldOperations<X,GenericTrait<X>,R>
//    , ProvideConcreteGenericComparisonOperations<X,GenericTrait<X>,LessTrait<X>,EqualsTrait<X>>
{
    friend OutputStream& operator<<(OutputStream& os, X const& x) { return Operations<X>::_write(os,x); }
    friend InputStream& operator>>(InputStream& is, X& x) { return Operations<X>::_read(is,x); }
};

template<class X> struct DispatchDirectedFloatOperations
    : DispatchDirectedNumericOperations<X,OppositeTrait<X>>
    , DispatchDirectedNumericOperations<OppositeTrait<X>,X>
    , DispatchDirectedComparisonOperations<X,OppositeTrait<X>,LessTrait<X>,EqualsTrait<X>>
    , DispatchDirectedComparisonOperations<OppositeTrait<X>,X,LessTrait<OppositeTrait<X>>,EqualsTrait<OppositeTrait<X>>>
    , DefineConcreteGenericOperators<X>
//    , ProvideConcreteGenericDirectedGroupOperations<X,OppositeTrait<X>,GenericTrait<X>,GenericTrait<OppositeTrait<X>>>
//    , ProvideConcreteGenericDirectedGroupOperations<OppositeTrait<X>,X,GenericTrait<OppositeTrait<X>>,GenericTrait<X>>
//    , ProvideConcreteGenericDirectedComparisonOperations<X,GenericTrait<OppositeTrait<X>>,LessTrait<X>,EqualsTrait<X>>
//    , ProvideConcreteGenericDirectedComparisonOperations<OppositeTrait<X>,GenericTrait<X>,LessTrait<OppositeTrait<X>>,EqualsTrait<OppositeTrait<X>>>
{
    friend OutputStream& operator<<(OutputStream& os, X const& x) { return Operations<X>::_write(os,x); }
    friend InputStream& operator>>(InputStream& is, X& x) { return Operations<X>::_read(is,x); }
};

template<class PX, class QPX=OppositeTrait<PX>> struct DispatchPositiveDirectedFloatOperations
    : public DispatchPositiveDirectedNumericOperations<PX,QPX>
    , public DispatchPositiveDirectedNumericOperations<QPX,PX>
    , DefineConcreteGenericOperators<PX>
//    , public ProvideConcreteGenericDirectedSemiFieldOperations<PX,QPX,Nat,Nat>
//    , public ProvideConcreteGenericDirectedSemiFieldOperations<QPX,PX,Nat,Nat>
{
};

template<class PX> struct DispatchPositiveFloatOperations
    : public DispatchPositiveDirectedNumericOperations<PX,PX>
    , DefineConcreteGenericOperators<PX>
//    , public ProvideConcreteGenericDirectedSemiFieldOperations<PX,PX,Nat,Nat>
{
};


template<class PR, class P1, class P2> using FloatWeakerType = FloatType<Weaker<P1,P2>,PR>;

template<class PR, class P> using NegatedFloatType = FloatType<Negated<P>,PR>;
template<class PR, class P> using FloatNegateType = FloatType<Negated<P>,PR>;

template<class PR, class P1, class P2> using FloatSumType = FloatType<Widen<Weaker<P1,P2>>,PR>;
template<class PR, class P1, class P2> using FloatDifferenceType = FloatType<Widen<Weaker<P1,Negated<P2>>>,PR>;
template<class PR, class P1, class P2> using FloatProductType = FloatType<Widen<Weaker<P1,P2>>,PR>;
template<class PR, class P1, class P2> using FloatQuotientType = FloatType<Widen<Weaker<P1,Inverted<P2>>>,PR>;

template<class PR, class P1, class P2> using FloatEqualsType = LogicalType<Equality<Weaker<P1,Negated<P2>>>>;
template<class PR, class P1, class P2> using FloatLessType = LogicalType<Generic<Weaker<P1,Negated<P2>>>>;

} // namespace Ariadne

#include "float_factory.hpp"

#endif
