/***************************************************************************
 *            numeric/real.h
 *
 *  Copyright 2013-14  Pieter Collins
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

/*! \file numeric/real.h
 *  \brief
 */



#ifndef ARIADNE_REAL_H
#define ARIADNE_REAL_H

#include "utility/typedefs.h"
#include "utility/pointer.h"

#include "numeric/logical.decl.h"
#include "numeric/number.decl.h"
#include "numeric/float.decl.h"

#include "numeric/arithmetic.h"

namespace Ariadne {

class Real;
class LowerReal;
class UpperReal;
class PositiveReal;
class PositiveLowerReal;
class PositiveUpperReal;
template<> struct IsNumericType<Real> : True { };

struct Accuracy {
    Nat _bits;
    Accuracy(Nat bits) : _bits(bits) { }
    Nat bits() const { return _bits; }
    TwoExp error() const;
    friend OutputStream& operator<<(OutputStream& os, Accuracy acc) { return os << "Accuracy("<<acc.bits()<<")"; }
};

extern const Real pi;
extern const Real infinity;

//! \ingroup UserNumericTypeSubModule
//! \brief Computable real numbers definable in terms of elementary functions.
class Real
    : public DeclareRealOperations<Real,PositiveReal>
    , public DeclareAnalyticFieldOperations<Real>
    , public DeclareLatticeOperations<Real,PositiveReal>
    , public DeclareComparisonOperations<Real,Kleenean,NegatedSierpinskian>
    , public DefineFieldOperators<Real>
{
  private: public:
    class Interface;
  private: public:
    SharedPointer<Interface> _ptr;
  private: public:
    explicit Real(SharedPointer<Interface>);
    explicit Real(double,double,double);
  public:
    typedef EffectiveTag Paradigm;
    typedef Real NumericType;
  public:
    Real();

    explicit Real(double);

    template<class M, EnableIf<And<IsIntegral<M>,IsUnsigned<M>>> = dummy> Real(M m);
    template<class N, EnableIf<And<IsIntegral<N>,IsSigned<N>>> = dummy> Real(N n);

    Real(Integer const& n);
    Real(Dyadic const& d);
    Real(Decimal const& d);
    Real(Rational const& q);

    explicit Real(Float64Value x); //!< DEPRECATED
    explicit Real(EffectiveNumber r); //!< DEPRECATED

    operator Number<EffectiveTag>() const;

    // Extract bounds
    UpperReal upper() const;
    LowerReal lower() const;
    Float64Approximation approx() const;
    double get_d() const;

    // Extract arbitrarily accurate approximations
    Float64Bounds operator() (Precision64 pr) const;
    FloatMPBounds operator() (PrecisionMP pr) const;
    Float64Bounds get(Precision64 pr) const;
    FloatMPBounds get(PrecisionMP pr) const;
    FloatMPBounds evaluate(Accuracy acc) const;


    friend PositiveReal abs(Real const&);

    friend PositiveUpperReal mag(Real const&);
    friend Float64Error mag(Real const&, Precision64);

    friend PositiveReal dist(Real const&, Real const&);

    friend Bool same(Real const&, Real const&);

    friend OutputStream& operator<<(OutputStream&, Real const&);

    friend Number<EffectiveTag> operator+(Number<EffectiveTag>, Number<EffectiveTag>);
    friend Number<EffectiveTag> operator-(Number<EffectiveTag>, Number<EffectiveTag>);
    friend Number<EffectiveTag> operator*(Number<EffectiveTag>, Number<EffectiveTag>);
    friend Number<EffectiveTag> operator/(Number<EffectiveTag>, Number<EffectiveTag>);
  private:
    Real(std::int64_t n, Void*);
    Real(std::uint64_t m, Void*);
};

template<class M, EnableIf<And<IsIntegral<M>,IsUnsigned<M>>>> inline Real::Real(M m) : Real(std::uint64_t(m),nullptr) { };
template<class N, EnableIf<And<IsIntegral<N>,IsSigned<N>>>> inline Real::Real(N n) : Real(std::int64_t(n),nullptr) { };


//! \ingroup UserNumericTypeSubModule
//! \brief Computable lower real numbers defined by conversion to concrete floats.
class LowerReal
    : public DirectedRing<LowerReal,UpperReal,PositiveLowerReal>
{
  private: public:
    SharedPointer<Real::Interface> _ptr;
  private: public:
    explicit LowerReal(SharedPointer<Real::Interface>);
  public:
    typedef EffectiveLowerTag Paradigm;
    typedef LowerReal NumericType;
  public:
    LowerReal(Real);
  public:
    Float64LowerBound operator() (Precision64 pr) const;
    FloatMPLowerBound operator() (PrecisionMP pr) const;
    Float64LowerBound get(Precision64 pr) const;
    FloatMPLowerBound get(PrecisionMP pr) const;
  public:
    friend LowerReal max(LowerReal const&, LowerReal const&);
    friend LowerReal min(LowerReal const&, LowerReal const&);

    friend LowerReal neg(UpperReal const&);
    friend UpperReal neg(LowerReal const&);
    friend LowerReal add(LowerReal const&, LowerReal const&);
    friend LowerReal sub(LowerReal const&, UpperReal const&);
    friend UpperReal sub(UpperReal const&, LowerReal const&);

    friend PositiveLowerReal sqrt(PositiveLowerReal const&);
    friend PositiveLowerReal exp(LowerReal const&);
    friend LowerReal log(PositiveLowerReal const&);
    friend LowerReal atan(LowerReal const&);

    friend PositiveLowerReal add(PositiveLowerReal const&, PositiveLowerReal const&);
    friend PositiveLowerReal mul(PositiveLowerReal const&, PositiveLowerReal const&);
};

//! \ingroup UserNumericTypeSubModule
//! \brief Computable lower real numbers defined by conversion to concrete floats.
class UpperReal
    : public DirectedRing<UpperReal,LowerReal,PositiveUpperReal>
{
  private: public:
    SharedPointer<Real::Interface> _ptr;
  private: public:
    explicit UpperReal(SharedPointer<Real::Interface>);
  public:
    typedef EffectiveUpperTag Paradigm;
    typedef UpperReal NumericType;
  public:
    UpperReal(Real);
  public:
    Float64UpperBound operator() (Precision64 pr) const;
    FloatMPUpperBound operator() (PrecisionMP pr) const;
    Float64UpperBound get(Precision64 pr) const;
    FloatMPUpperBound get(PrecisionMP pr) const;
  public:
    friend UpperReal max(UpperReal const&, UpperReal const&);
    friend UpperReal min(UpperReal const&, UpperReal const&);

    friend UpperReal neg(LowerReal const&);
    friend LowerReal neg(UpperReal const&);
    friend UpperReal add(UpperReal const&, UpperReal const&);
    friend UpperReal sub(UpperReal const&, LowerReal const&);
    friend LowerReal sub(LowerReal const&, UpperReal const&);

    friend PositiveUpperReal sqrt(PositiveUpperReal const&);
    friend PositiveUpperReal exp(UpperReal const&);
    friend UpperReal log(PositiveUpperReal const&);
    friend UpperReal atan(UpperReal const&);

    friend PositiveUpperReal add(PositiveUpperReal const&, PositiveUpperReal const&);
    friend PositiveUpperReal mul(PositiveUpperReal const&, PositiveUpperReal const&);
};

//! \ingroup UserNumericTypeSubModule
//! \brief Computable lower real numbers defined by conversion to concrete floats.
class PositiveReal : public Real
{
  public:
    using Real::Real;
    PositiveReal() : Real() { }
    PositiveReal(Real r) : Real(r) { }
    PositiveFloat64Bounds get(Precision64 pr) const;
    PositiveFloatMPBounds get(PrecisionMP pr) const;
  public:
    PositiveReal max(PositiveReal const&, PositiveReal const&);
    PositiveReal max(Real const&, PositiveReal const&);
    PositiveReal max(PositiveReal const&, Real const&);

    PositiveReal min(PositiveReal const&, PositiveReal const&);

    PositiveReal rec(PositiveReal const&);
    PositiveReal add(PositiveReal const&, PositiveReal const&);
    PositiveReal mul(PositiveReal const&, PositiveReal const&);
    PositiveReal div(PositiveReal const&, PositiveReal const&);
};

//! \ingroup UserNumericTypeSubModule
//! \brief Computable lower real numbers defined by conversion to concrete floats.
class PositiveLowerReal : public LowerReal
{
  public:
    using LowerReal::LowerReal;
    PositiveLowerReal(LowerReal r) : LowerReal(r) { }
    PositiveFloat64LowerBound get(Precision64 pr) const;
    PositiveFloatMPLowerBound get(PrecisionMP pr) const;
  public:
    PositiveLowerReal rec(PositiveUpperReal const&);
    PositiveUpperReal rec(PositiveLowerReal const&);
    PositiveLowerReal add(PositiveLowerReal const&, PositiveLowerReal const&);
    PositiveLowerReal mul(PositiveLowerReal const&, PositiveLowerReal const&);
    PositiveLowerReal div(PositiveLowerReal const&, PositiveUpperReal const&);
    PositiveUpperReal div(PositiveUpperReal const&, PositiveLowerReal const&);
};

//! \ingroup UserNumericTypeSubModule
//! \brief Computable lower real numbers defined by conversion to concrete floats.
class PositiveUpperReal : public UpperReal
{
  public:
    using UpperReal::UpperReal;
    PositiveUpperReal(UpperReal r) : UpperReal(r) { }
    PositiveFloat64UpperBound get(Precision64 pr) const;
    PositiveFloatMPUpperBound get(PrecisionMP pr) const;
  public:
    PositiveUpperReal rec(PositiveLowerReal const&);
    PositiveLowerReal rec(PositiveUpperReal const&);
    PositiveUpperReal add(PositiveUpperReal const&, PositiveUpperReal const&);
    PositiveUpperReal mul(PositiveUpperReal const&, PositiveUpperReal const&);
    PositiveUpperReal div(PositiveUpperReal const&, PositiveLowerReal const&);
    PositiveLowerReal div(PositiveLowerReal const&, PositiveUpperReal const&);
};

PositiveReal cast_positive(Real const& x);

} // namespace Ariadne

#include "numeric/logical.h"

namespace Ariadne {

template<class N, EnableIf<IsIntegral<N>> =dummy> inline Kleenean operator==(const Real& x1, N n2) { return x1==Real(n2); }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline Kleenean operator!=(const Real& x1, N n2) { return x1!=Real(n2); }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline Kleenean operator<=(const Real& x1, N n2) { return x1<=Real(n2); }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline Kleenean operator>=(const Real& x1, N n2) { return x1>=Real(n2); }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline Kleenean operator< (const Real& x1, N n2) { return x1< Real(n2); }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline Kleenean operator> (const Real& x1, N n2) { return x1> Real(n2); }

/*
template<class D, EnableIf<IsFloatingPoint<D>> =dummy> inline auto
    operator==(Real r, D d) -> decltype(r==Float64Approximation(d)) { return r==Float64Approximation(d); }
template<class D, EnableIf<IsFloatingPoint<D>> =dummy> inline auto
    operator!=(Real r, D d) -> decltype(r!=Float64Approximation(d)) { return r!=Float64Approximation(d); }
template<class D, EnableIf<IsFloatingPoint<D>> =dummy> inline auto
    operator< (Real r, D d) -> decltype(r< Float64Approximation(d)) { return r< Float64Approximation(d); }
template<class D, EnableIf<IsFloatingPoint<D>> =dummy> inline auto
    operator> (Real r, D d) -> decltype(r> Float64Approximation(d)) { return r> Float64Approximation(d); }
template<class D, EnableIf<IsFloatingPoint<D>> =dummy> inline auto
    operator<=(Real r, D d) -> decltype(r<=Float64Approximation(d)) { return r<=Float64Approximation(d); }
template<class D, EnableIf<IsFloatingPoint<D>> =dummy> inline auto
    operator>=(Real r, D d) -> decltype(r>=Float64Approximation(d)) { return r>=Float64Approximation(d); }
*/
/*
template<class T> auto operator+(T const& t) -> decltype(pos(t)) { return pos(t); }
template<class T> auto operator-(T const& t) -> decltype(neg(t)) { return neg(t); }
template<class T1, class T2> auto operator+(T1 const& t1, T2 const& t2) -> decltype(add(t1,t2)) { return add(t1,t2); }
template<class T1, class T2> auto operator-(T1 const& t1, T2 const& t2) -> decltype(sub(t1,t2)) { return sub(t1,t2); }
template<class T1, class T2> auto operator*(T1 const& t1, T2 const& t2) -> decltype(mul(t1,t2)) { return mul(t1,t2); }
template<class T1, class T2> auto operator/(T1 const& t1, T2 const& t2) -> decltype(div(t1,t2)) { return div(t1,t2); }
*/
} // namespace Ariadne

#endif
