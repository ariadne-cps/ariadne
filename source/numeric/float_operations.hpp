/***************************************************************************
 *            numeric/float_operations.hpp
 *
 *  Copyright  2008-20  Pieter Collins
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

/*! \file numeric/float_operations.hpp
 *  \brief Common header for floating-point number operations.
 */

#ifndef ARIADNE_FLOAT_OPERATIONS_HPP
#define ARIADNE_FLOAT_OPERATIONS_HPP

#include "utility/macros.hpp"

#include "number.decl.hpp"
#include "float.decl.hpp"
#include "arithmetic.hpp"

#include "sign.hpp"
#include "operators.hpp"
#include "logical.hpp"

namespace Ariadne {

template<class X> X generic_pow(X const& x, Nat m) {
    X r=x; r*=0; r+=1; X p=x; while(m!=0) { if(m%2==1) { r=r*p; } p=p*p; m=m/2; } return r; }
template<class X> X generic_pow(const X& x, Int n) {
    return n>=0 ? generic_pow(x,Nat(n)) : rec(generic_pow(x,Nat(-n))); }


template<class... PRS> class Float;
template<class F, class FE> class Ball;
template<class F> class Bounds;
template<class F> class UpperBound;
template<class F> class LowerBound;
template<class F> class Approximation;
template<class F> class Error;
template<class F> class Rounded;

template<class X> struct IsConcreteNumber : False { };
template<class PR> struct IsConcreteNumber<Float<PR>> : True { };
template<class PR> struct IsConcreteNumber<Rounded<Float<PR>>> : True { };
template<class PR,class PRE> struct IsConcreteNumber<Ball<Float<PR>,Float<PRE>>> : True { };
template<class PR> struct IsConcreteNumber<Bounds<Float<PR>>> : True { };
template<class PR> struct IsConcreteNumber<UpperBound<Float<PR>>> : True { };
template<class PR> struct IsConcreteNumber<LowerBound<Float<PR>>> : True { };
template<class PR> struct IsConcreteNumber<Approximation<Float<PR>>> : True { };
template<class PR> struct IsConcreteNumber<Error<Float<PR>>> : True { };
template<class X> struct IsConcreteNumber<Positive<X>> : IsConcreteNumber<X> { };

template<class X> concept AConcreteNumber = IsConcreteNumber<X>::value;

template<AConcreteNumber X, GenericNumber Y> inline decltype(auto) add(X const& x, Y const& y) { return add(x,factory(x).create(y)); }
template<AConcreteNumber X, GenericNumber Y> inline decltype(auto) sub(X const& x, Y const& y) { return sub(x,factory(x).create(y)); }
template<AConcreteNumber X, GenericNumber Y> inline decltype(auto) mul(X const& x, Y const& y) { return mul(x,factory(x).create(y)); }
template<AConcreteNumber X, GenericNumber Y> inline decltype(auto) div(X const& x, Y const& y) { return div(x,factory(x).create(y)); }
template<AConcreteNumber X, GenericNumber Y> inline decltype(auto) add(Y const& y, X const& x) { return add(factory(x).create(y),x); }
template<AConcreteNumber X, GenericNumber Y> inline decltype(auto) sub(Y const& y, X const& x) { return sub(factory(x).create(y),x); }
template<AConcreteNumber X, GenericNumber Y> inline decltype(auto) mul(Y const& y, X const& x) { return mul(factory(x).create(y),x); }
template<AConcreteNumber X, GenericNumber Y> inline decltype(auto) div(Y const& y, X const& x) { return div(factory(x).create(y),x); }

template<AConcreteNumber X, GenericNumber Y> inline decltype(auto) max(X const& x, Y const& y) {
    if constexpr (Convertible<Y,Dyadic>) { return max(x,factory(x).create(y,ExactTag())); } else { return max(x,factory(x).create(y)); } }
template<AConcreteNumber X, GenericNumber Y> inline decltype(auto) min(X const& x, Y const& y) {
    if constexpr (Convertible<Y,Dyadic>) { return min(x,factory(x).create(y,ExactTag())); } else { return min(x,factory(x).create(y)); } }
template<AConcreteNumber X, GenericNumber Y> inline decltype(auto) max(Y const& y, X const& x) {
    if constexpr (Convertible<Y,Dyadic>) { return max(factory(x).create(y,ExactTag()),x); } else { return max(factory(x).create(y),x); } }
template<AConcreteNumber X, GenericNumber Y> inline decltype(auto) min(Y const& y, X const& x) {
    if constexpr (Convertible<Y,Dyadic>) { return min(factory(x).create(y,ExactTag()),x); } else { return min(factory(x).create(y),x); } }

template<AConcreteNumber X, GenericNumber Y> inline decltype(auto) operator+(X const& x, Y const& y) { return x+factory(x).create(y); }
template<AConcreteNumber X, GenericNumber Y> inline decltype(auto) operator-(X const& x, Y const& y) { return x-factory(x).create(y); }
template<AConcreteNumber X, GenericNumber Y> inline decltype(auto) operator*(X const& x, Y const& y) { return x*factory(x).create(y); }
template<AConcreteNumber X, GenericNumber Y> inline decltype(auto) operator/(X const& x, Y const& y) { return x/factory(x).create(y); }

template<AConcreteNumber X, GenericNumber Y> inline decltype(auto) operator+(Y const& y, X const& x) { return factory(x).create(y)+x; }
template<AConcreteNumber X, GenericNumber Y> inline decltype(auto) operator-(Y const& y, X const& x) { return factory(x).create(y)-x; }
template<AConcreteNumber X, GenericNumber Y> inline decltype(auto) operator*(Y const& y, X const& x) { return factory(x).create(y)*x; }
template<AConcreteNumber X, GenericNumber Y> inline decltype(auto) operator/(Y const& y, X const& x) { return factory(x).create(y)/x; }

template<AConcreteNumber X, GenericNumber Y> inline decltype(auto) operator+=(X& x, Y const& y) { return x+=factory(x).create(y); }
template<AConcreteNumber X, GenericNumber Y> inline decltype(auto) operator-=(X& x, Y const& y) { return x-=factory(x).create(y); }
template<AConcreteNumber X, GenericNumber Y> inline decltype(auto) operator*=(X& x, Y const& y) { return x*=factory(x).create(y); }
template<AConcreteNumber X, GenericNumber Y> inline decltype(auto) operator/=(X& x, Y const& y) { return x/=factory(x).create(y); }

template<class Y1, class Y2> concept HaveCmp = requires(Y1 y1, Y2 y2) {
    { cmp(y1,y2) } -> SameAs<Comparison>;
};

/*
template<class OP, class Y1, class Y2> requires HaveCmp<Y1,Y2> Boolean _cmp(OP op, Y1 const& y1, Y2 const& y2) {
    return op(cmp(y1,y2),Comparison::EQUAL); }
template<class Y1, class Y2> requires HaveCmp<Y1,Y2> Boolean operator==(Y1 const& y1, Y2 const& y2) { return _cmp(Eq(),y1,y2); }
template<class Y1, class Y2> requires HaveCmp<Y1,Y2> Boolean operator!=(Y1 const& y1, Y2 const& y2) { return _cmp(Ne(),y1,y2); }
template<class Y1, class Y2> requires HaveCmp<Y1,Y2> Boolean operator< (Y1 const& y1, Y2 const& y2) { return _cmp(Lt(),y1,y2); }
template<class Y1, class Y2> requires HaveCmp<Y1,Y2> Boolean operator> (Y1 const& y1, Y2 const& y2) { return _cmp(Gt(),y1,y2); }
template<class Y1, class Y2> requires HaveCmp<Y1,Y2> Boolean operator<=(Y1 const& y1, Y2 const& y2) { return _cmp(Le(),y1,y2); }
template<class Y1, class Y2> requires HaveCmp<Y1,Y2> Boolean operator>=(Y1 const& y1, Y2 const& y2) { return _cmp(Ge(),y1,y2); }
*/

template<class Y1, class Y2> requires HaveCmp<Y1,Y2> Boolean operator==(Y1 const& y1, Y2 const& y2) { return cmp(y1,y2)==Comparison::EQUAL; }
template<class Y1, class Y2> requires HaveCmp<Y1,Y2> Boolean operator!=(Y1 const& y1, Y2 const& y2) { return cmp(y1,y2)!=Comparison::EQUAL; }
template<class Y1, class Y2> requires HaveCmp<Y1,Y2> Boolean operator< (Y1 const& y1, Y2 const& y2) { return cmp(y1,y2)< Comparison::EQUAL; }
template<class Y1, class Y2> requires HaveCmp<Y1,Y2> Boolean operator> (Y1 const& y1, Y2 const& y2) { return cmp(y1,y2)> Comparison::EQUAL; }
template<class Y1, class Y2> requires HaveCmp<Y1,Y2> Boolean operator<=(Y1 const& y1, Y2 const& y2) { return cmp(y1,y2)<=Comparison::EQUAL; }
template<class Y1, class Y2> requires HaveCmp<Y1,Y2> Boolean operator>=(Y1 const& y1, Y2 const& y2) { return cmp(y1,y2)>=Comparison::EQUAL; }

template<AConcreteNumber X, GenericNumber Y> requires (not HaveCmp<X,Y>) inline decltype(auto) operator==(X const& x, Y const& y) {
    return x==factory(x).create(y); }
template<AConcreteNumber X, GenericNumber Y> requires (not HaveCmp<X,Y>) inline decltype(auto) operator!=(X const& x, Y const& y) {
    return x!=factory(x).create(y); }
template<AConcreteNumber X, GenericNumber Y> requires (not HaveCmp<X,Y>) inline decltype(auto) operator<=(X const& x, Y const& y) {
     return x<=factory(x).create(y); }
template<AConcreteNumber X, GenericNumber Y> requires (not HaveCmp<X,Y>) inline decltype(auto) operator>=(X const& x, Y const& y) {
     return x>=factory(x).create(y); }
template<AConcreteNumber X, GenericNumber Y> requires (not HaveCmp<X,Y>) inline decltype(auto) operator< (X const& x, Y const& y) {
     return x< factory(x).create(y); }
template<AConcreteNumber X, GenericNumber Y> requires (not HaveCmp<X,Y>) inline decltype(auto) operator> (X const& x, Y const& y) {
     return x> factory(x).create(y); }

template<AConcreteNumber X, GenericNumber Y> requires (not HaveCmp<X,Y>) inline decltype(auto) operator==(Y const& y, X const& x) {
     return factory(x).create(y)==x; }
template<AConcreteNumber X, GenericNumber Y> requires (not HaveCmp<X,Y>) inline decltype(auto) operator!=(Y const& y, X const& x) {
     return factory(x).create(y)!=x; }
template<AConcreteNumber X, GenericNumber Y> requires (not HaveCmp<X,Y>) inline decltype(auto) operator<=(Y const& y, X const& x) {
     return factory(x).create(y)<=x; }
template<AConcreteNumber X, GenericNumber Y> requires (not HaveCmp<X,Y>) inline decltype(auto) operator>=(Y const& y, X const& x) {
     return factory(x).create(y)>=x; }
template<AConcreteNumber X, GenericNumber Y> requires (not HaveCmp<X,Y>) inline decltype(auto) operator< (Y const& y, X const& x) {
     return factory(x).create(y)< x; }
template<AConcreteNumber X, GenericNumber Y> requires (not HaveCmp<X,Y>) inline decltype(auto) operator> (Y const& y, X const& x) {
     return factory(x).create(y)> x; }


template<class X, class R=X> struct DeclareRealOperations { };
template<class X, class LT, class EQ> struct DeclareComparisonOperations { };
template<class X, class Y, class R> struct DeclareMixedFieldOperators { };
template<class X, class NX> struct DeclareDirectedNumericOperations { };
template<class X, class NX, class Y, class NY> struct DeclareMixedDirectedGroupOperators { };
template<class PX, class QPX> struct DeclarePositiveDirectedNumericOperations { };
template<class X, class R=X> struct DispatchNumericOperations { };
template<class X, class NX, class R=X> struct DispatchDirectedNumericOperations { };
template<class PX, class QPX> struct DispatchPositiveDirectedNumericOperations { };
template<class X, class LT, class EQ> struct DispatchComparisonOperations { };
template<class X, class NX, class LT, class EQ> struct DispatchDirectedComparisonOperations { };
template<class X, class R=X> struct DefineFieldOperators { };
template<class X, class NX> class DefineDirectedGroupOperators { };
template<class X, class LT, class EQ=LT> struct DefineComparisonOperators { };
template<class X, class NX, class LT, class EQ> struct DefineDirectedComparisonOperators { };
template<class X> struct DefineConcreteGenericOperators { };

template<class X, class Y=Void, class R=X> struct ProvideConcreteGenericElementaryOperations { };





template<class X, class Y=Void, class R=X> struct ProvideConcreteGenericFieldOperators { };
template<class X, class Y=Void, class R=X> struct ProvideConcreteGenericFieldOperations { };

template<class X> struct ProvideConcreteGenericFieldOperators<X,Void> {
    template<GenericNumber Y> friend decltype(auto) operator+(X const& x, Y const& y) { return add(x,factory(x).create(y)); }
    template<GenericNumber Y> friend decltype(auto) operator-(X const& x, Y const& y) { return sub(x,factory(x).create(y)); }
    template<GenericNumber Y> friend decltype(auto) operator*(X const& x, Y const& y) { return mul(x,factory(x).create(y)); }
    template<GenericNumber Y> friend decltype(auto) operator/(X const& x, Y const& y) { return div(x,factory(x).create(y)); }
    template<GenericNumber Y> friend decltype(auto) operator+(Y const& y, X const& x) { return add(factory(x).create(y),x); }
    template<GenericNumber Y> friend decltype(auto) operator-(Y const& y, X const& x) { return sub(factory(x).create(y),x); }
    template<GenericNumber Y> friend decltype(auto) operator*(Y const& y, X const& x) { return mul(factory(x).create(y),x); }
    template<GenericNumber Y> friend decltype(auto) operator/(Y const& y, X const& x) { return div(factory(x).create(y),x); }
    template<GenericNumber Y> friend decltype(auto) operator+=(X& x, Y const& y) { return x=add(x,factory(x).create(y)); }
    template<GenericNumber Y> friend decltype(auto) operator-=(X& x, Y const& y) { return x=sub(x,factory(x).create(y)); }
    template<GenericNumber Y> friend decltype(auto) operator*=(X& x, Y const& y) { return x=mul(x,factory(x).create(y)); }
    template<GenericNumber Y> friend decltype(auto) operator/=(X& x, Y const& y) { return x=div(x,factory(x).create(y)); }
};
template<class X> struct ProvideConcreteGenericFieldOperations<X,Void> : ProvideConcreteGenericFieldOperators<X,Void> {
    template<GenericNumber Y> friend decltype(auto) add(X const& x, Y const& y) { return add(x,factory(x).create(y)); }
    template<GenericNumber Y> friend decltype(auto) sub(X const& x, Y const& y) { return sub(x,factory(x).create(y)); }
    template<GenericNumber Y> friend decltype(auto) mul(X const& x, Y const& y) { return mul(x,factory(x).create(y)); }
    template<GenericNumber Y> friend decltype(auto) div(X const& x, Y const& y) { return div(x,factory(x).create(y)); }
    template<GenericNumber Y> friend decltype(auto) add(Y const& y, X const& x) { return add(factory(x).create(y),x); }
    template<GenericNumber Y> friend decltype(auto) sub(Y const& y, X const& x) { return sub(factory(x).create(y),x); }
    template<GenericNumber Y> friend decltype(auto) mul(Y const& y, X const& x) { return mul(factory(x).create(y),x); }
    template<GenericNumber Y> friend decltype(auto) div(Y const& y, X const& x) { return div(factory(x).create(y),x); }
};

template<class X, class Y=Void, class R=X> struct ProvideConcreteGenericArithmeticOperators : ProvideConcreteGenericFieldOperators<X,Y,R> { };
template<class X, class Y=Void, class R=X> struct ProvideConcreteGenericArithmeticOperations : ProvideConcreteGenericFieldOperations<X,Y,R> { };


template<class X> class Operations;

template<class X, class R=X> struct DeclareFloatOperations
    : DeclareRealOperations<X,R>
    , DeclareComparisonOperations<X,LessTrait<X>,EqualsTrait<X>>
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
{
};

template<class PX> struct DeclarePositiveFloatOperations
    : DeclarePositiveDirectedFloatOperations<PX,PX>
{
};

template<class X, class R=X> struct DefineFloatOperations
    : public DefineFieldOperators<X,R>
    , public DefineComparisonOperators<X,LessTrait<X>,EqualsTrait<X>>
    , public DefineConcreteGenericOperators<X>
    , public DeclareFloatOperations<X>
{
};

template<class X, class NX> struct DefineDirectedFloatOperations
    : public DefineDirectedGroupOperators<X,NX>
    , public DefineDirectedGroupOperators<NX,X>
    , public DefineDirectedComparisonOperators<X,NX,LessTrait<X>,EqualsTrait<X>>
    , public DefineDirectedComparisonOperators<NX,X,LessTrait<NX>,EqualsTrait<NX>>
    , public DefineConcreteGenericOperators<X>
{
};

template<class X, class R=X> struct DispatchFloatOperations
    : DispatchNumericOperations<X,R>
    , DispatchComparisonOperations<X,LessTrait<X>,EqualsTrait<X>>
    , DefineConcreteGenericOperators<X>
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
{
    friend OutputStream& operator<<(OutputStream& os, X const& x) { return Operations<X>::_write(os,x); }
    friend InputStream& operator>>(InputStream& is, X& x) { return Operations<X>::_read(is,x); }
};

template<class PX, class QPX=OppositeTrait<PX>> struct DispatchPositiveDirectedFloatOperations
    : DispatchPositiveDirectedNumericOperations<PX,QPX>
    , DispatchPositiveDirectedNumericOperations<QPX,PX>
    , DefineConcreteGenericOperators<PX>
{
};

template<class PX> struct DispatchPositiveFloatOperations
    : DispatchPositiveDirectedNumericOperations<PX,PX>
    , DefineConcreteGenericOperators<PX>
{
};


/*
// Mixed Float-ValidatedNumber comparisons
template<ARawFloat F, ANonExactValidatedNumber Y> ValidatedKleenean operator==(Y const& y1, F const& x2) { return Bounds<F>(y1,x2.precision())==x2; }
template<ARawFloat F, ANonExactValidatedNumber Y> ValidatedKleenean operator!=(Y const& y1, F const& x2) { return Bounds<F>(y1,x2.precision())!=x2; }
template<ARawFloat F, ANonExactValidatedNumber Y> ValidatedKleenean operator< (Y const& y1, F const& x2) { return Bounds<F>(y1,x2.precision())< x2; }
template<ARawFloat F, ANonExactValidatedNumber Y> ValidatedKleenean operator> (Y const& y1, F const& x2) { return Bounds<F>(y1,x2.precision())> x2; }
template<ARawFloat F, ANonExactValidatedNumber Y> ValidatedKleenean operator<=(Y const& y1, F const& x2) { return Bounds<F>(y1,x2.precision())<=x2; }
template<ARawFloat F, ANonExactValidatedNumber Y> ValidatedKleenean operator>=(Y const& y1, F const& x2) { return Bounds<F>(y1,x2.precision())>=x2; }
template<ARawFloat F, ANonExactValidatedNumber Y> ValidatedKleenean operator==(F const& x1, Y const& y2) { return x1==Bounds<F>(y2,x1.precision()); }
template<ARawFloat F, ANonExactValidatedNumber Y> ValidatedKleenean operator!=(F const& x1, Y const& y2) { return x1!=Bounds<F>(y2,x1.precision()); }
template<ARawFloat F, ANonExactValidatedNumber Y> ValidatedKleenean operator< (F const& x1, Y const& y2) { return x1< Bounds<F>(y2,x1.precision()); }
template<ARawFloat F, ANonExactValidatedNumber Y> ValidatedKleenean operator> (F const& x1, Y const& y2) { return x1> Bounds<F>(y2,x1.precision()); }
template<ARawFloat F, ANonExactValidatedNumber Y> ValidatedKleenean operator<=(F const& x1, Y const& y2) { return x1<=Bounds<F>(y2,x1.precision()); }
template<ARawFloat F, ANonExactValidatedNumber Y> ValidatedKleenean operator>=(F const& x1, Y const& y2) { return x1>=Bounds<F>(y2,x1.precision()); }

// Mixed Float-ValidatedNumber operations
template<ARawFloat F, AValidatedNumber Y> Bounds<F> operator+(F const& x1, Y const& y2) {
    if constexpr (AnExactDyadic<Y>) { return operator+(x1,F(y2,x1.precision())); }
    else { return operator+(x1,Bounds<F>(y2,x1.precision())); } }
template<ARawFloat F, AValidatedNumber Y> Bounds<F> operator-(F const& x1, Y const& y2) {
    if constexpr (AnExactDyadic<Y>) { return operator-(x1,F(y2,x1.precision())); }
    else { return operator-(x1,Bounds<F>(y2,x1.precision())); } }
template<ARawFloat F, AValidatedNumber Y> Bounds<F> operator*(F const& x1, Y const& y2) {
    if constexpr (AnExactDyadic<Y>) { return operator*(x1,F(y2,x1.precision())); }
    else { return operator*(x1,Bounds<F>(y2,x1.precision())); } }
template<ARawFloat F, AValidatedNumber Y> Bounds<F> operator/(F const& x1, Y const& y2) {
    if constexpr (AnExactDyadic<Y>) { return operator/(x1,F(y2,x1.precision())); }
    else { return operator/(x1,Bounds<F>(y2,x1.precision())); } }
template<ARawFloat F, AValidatedNumber Y> Bounds<F> operator+(Y const& y1, F const& x2) {
    if constexpr (AnExactDyadic<Y>) { return operator+(F(y1,x2.precision()),x2); }
    else { return operator+(Bounds<F>(y1,x2.precision()),x2); } }
template<ARawFloat F, AValidatedNumber Y> Bounds<F> operator-(Y const& y1, F const& x2) {
    if constexpr (AnExactDyadic<Y>) { return operator-(F(y1,x2.precision()),x2); }
    else { return operator-(Bounds<F>(y1,x2.precision()),x2); } }
template<ARawFloat F, AValidatedNumber Y> Bounds<F> operator*(Y const& y1, F const& x2) {
    if constexpr (AnExactDyadic<Y>) { return operator*(F(y1,x2.precision()),x2); }
    else { return operator*(Bounds<F>(y1,x2.precision()),x2); } }
template<ARawFloat F, AValidatedNumber Y> Bounds<F> operator/(Y const& y1, F const& x2) {
    if constexpr (AnExactDyadic<Y>) { return operator/(F(y1,x2.precision()),x2); }
    else { return operator/(Bounds<F>(y1,x2.precision()),x2); } }

template<ARawFloat F, AValidatedNumber Y> Bounds<F> add(F const& x1, Y const& y2) {
    if constexpr (AnExactDyadic<Y>) { return add(x1,F(y2,x1.precision())); }
    else { return add(x1,Bounds<F>(y2,x1.precision())); } }
template<ARawFloat F, AValidatedNumber Y> Bounds<F> sub(F const& x1, Y const& y2) {
    if constexpr (AnExactDyadic<Y>) { return sub(x1,F(y2,x1.precision())); }
    else { return sub(x1,Bounds<F>(y2,x1.precision())); } }
template<ARawFloat F, AValidatedNumber Y> Bounds<F> mul(F const& x1, Y const& y2) {
    if constexpr (AnExactDyadic<Y>) { return mul(x1,F(y2,x1.precision())); }
    else { return mul(x1,Bounds<F>(y2,x1.precision())); } }
template<ARawFloat F, AValidatedNumber Y> Bounds<F> div(F const& x1, Y const& y2) {
    if constexpr (AnExactDyadic<Y>) { return div(x1,F(y2,x1.precision())); }
    else { return div(x1,Bounds<F>(y2,x1.precision())); } }
template<ARawFloat F, AValidatedNumber Y> Bounds<F> add(Y const& y1, F const& x2) {
    if constexpr (AnExactDyadic<Y>) { return add(F(y1,x2.precision()),x2); }
    else { return add(Bounds<F>(y1,x2.precision()),x2); } }
template<ARawFloat F, AValidatedNumber Y> Bounds<F> sub(Y const& y1, F const& x2) {
    if constexpr (AnExactDyadic<Y>) { return sub(F(y1,x2.precision()),x2); }
    else { return sub(Bounds<F>(y1,x2.precision()),x2); } }
template<ARawFloat F, AValidatedNumber Y> Bounds<F> mul(Y const& y1, F const& x2) {
    if constexpr (AnExactDyadic<Y>) { return mul(F(y1,x2.precision()),x2); }
    else { return mul(Bounds<F>(y1,x2.precision()),x2); } }
template<ARawFloat F, AValidatedNumber Y> Bounds<F> div(Y const& y1, F const& x2) {
    if constexpr (AnExactDyadic<Y>) { return div(F(y1,x2.precision()),x2); }
    else { return div(Bounds<F>(y1,x2.precision()),x2); } }

template<ARawFloat F, AValidatedNumber Y> decltype(auto) max(F const& x1, Y const& y2) {
    if constexpr (AnExactDyadic<Y>) { return max(x1,F(y2,x1.precision())); }
    else { return max(x1,Bounds<F>(y2,x1.precision())); } }
template<ARawFloat F, AValidatedNumber Y> decltype(auto) min(F const& x1, Y const& y2) {
    if constexpr (AnExactDyadic<Y>) { return min(x1,F(y2,x1.precision())); }
    else { return min(x1,Bounds<F>(y2,x1.precision())); } }
template<ARawFloat F, AValidatedNumber Y> decltype(auto) max(Y const& y1, F const& x2) {
    if constexpr (AnExactDyadic<Y>) { return max(F(y1,x2.precision()),x2); }
    else { return max(Bounds<F>(y1,x2.precision()),x2); } }
template<ARawFloat F, AValidatedNumber Y> decltype(auto) min(Y const& y1, F const& x2) {
    if constexpr (AnExactDyadic<Y>) { return min(F(y1,x2.precision()),x2); }
    else { return min(Bounds<F>(y1,x2.precision()),x2); } }
*/


template<class X> class Bounds;

template<class FLT> class DeclareRoundedOperations {
    friend Bounds<FloatDP> operator+(FloatDP const& x1, FloatDP const& x2);
    friend Bounds<FloatDP> operator-(FloatDP const& x1, FloatDP const& x2);
    friend Bounds<FloatDP> operator*(FloatDP const& x1, FloatDP const& x2);
    friend Bounds<FloatDP> operator/(FloatDP const& x1, FloatDP const& x2);

    friend Bounds<FloatDP> rec(FloatDP const& x);
    friend Bounds<FloatDP> add(FloatDP const& x1, FloatDP const& x2);
    friend Bounds<FloatDP> sub(FloatDP const& x1, FloatDP const& x2);
    friend Bounds<FloatDP> mul(FloatDP const& x1, FloatDP const& x2);
    friend Bounds<FloatDP> div(FloatDP const& x1, FloatDP const& x2);
    friend Bounds<FloatDP> fma(FloatDP const& x1, FloatDP const& x2, FloatDP const& x3);
    friend Bounds<FloatDP> pow(FloatDP const& x, Nat m);
    friend Bounds<FloatDP> pow(FloatDP const& x, Int n);
    friend Bounds<FloatDP> sqrt(FloatDP const& x);
    friend Bounds<FloatDP> exp(FloatDP const& x);
    friend Bounds<FloatDP> log(FloatDP const& x);
    friend Bounds<FloatDP> sin(FloatDP const& x);
    friend Bounds<FloatDP> cos(FloatDP const& x);
    friend Bounds<FloatDP> tan(FloatDP const& x);
    friend Bounds<FloatDP> asin(FloatDP const& x);
    friend Bounds<FloatDP> acos(FloatDP const& x);
    friend Bounds<FloatDP> atan(FloatDP const& x);
};

} // namespace Ariadne

#include "floatdp.hpp"
#include "floatmp.hpp"

#include "float_factory.hpp"

#endif
