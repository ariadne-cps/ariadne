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

    explicit Real(Float64Value x);

    operator Number<EffectiveTag>() const;

    // Extract floating-point properties
    Float64UpperBound upper() const;
    Float64LowerBound lower() const;
    Float64Approximation approx() const;
    double get_d() const;

    // Extract arbitrarily accurate approximations
    Float64Bounds operator() (Precision64 pr) const;
    FloatMPBounds operator() (PrecisionMP pr) const;
    Float64Bounds get(Precision64 pr) const;
    FloatMPBounds get(PrecisionMP pr) const;
    FloatMPBounds evaluate(Accuracy acc) const;

    // Non-templated to allow conversions
    friend Real add(Real,Real);
    friend Real sub(Real,Real);
    friend Real mul(Real,Real);
    friend Real div(Real,Real);
    friend Real pow(Real, Nat);
    friend Real pow(Real, Int);
    friend Real pos(Real);
    friend Real sqr(Real);
    friend Real neg(Real);
    friend Real rec(Real);
    friend Real sqrt(Real);
    friend Real exp(Real);
    friend Real log(Real);
    friend Real sin(Real);
    friend Real cos(Real);
    friend Real tan(Real);
    friend Real atan(Real);

    friend Real max(Real,Real);
    friend Real min(Real,Real);
    friend PositiveReal abs(Real);

    friend Float64Error mag(Real);

    friend ValidatedNegatedSierpinskian eq(Real,Real);
    friend ValidatedKleenean lt(Real,Real);

    friend PositiveReal dist(Real,Real);

    friend Real& operator+=(Real&,Real);
    friend Real& operator-=(Real&,Real);
    friend Real& operator*=(Real&,Real);
    friend Real& operator/=(Real&,Real);

    friend Real operator+(Real);
    friend Real operator-(Real);
    friend Real operator+(Real,Real);
    friend Real operator-(Real,Real);
    friend Real operator*(Real,Real);
    friend Real operator/(Real,Real);

    friend Falsifyable operator==(Real,Real);
    friend Verifyable operator!=(Real,Real);
    friend Quasidecidable operator<=(Real,Real);
    friend Quasidecidable operator>=(Real,Real);
    friend Quasidecidable operator< (Real,Real);
    friend Quasidecidable operator> (Real,Real);

    friend Bool same(Real,Real);

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
    LowerReal max(LowerReal,LowerReal);

    LowerReal min(LowerReal,LowerReal);
    Real min(Real,LowerReal);
    Real min(LowerReal,Real);

    LowerReal neg(UpperReal);
    UpperReal neg(LowerReal);
    LowerReal add(LowerReal,LowerReal);
    LowerReal sub(LowerReal,UpperReal);
    UpperReal sub(UpperReal,LowerReal);
};

//! \ingroup UserNumericTypeSubModule
//! \brief Computable lower real numbers defined by conversion to concrete floats.
class UpperReal
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
    UpperReal max(UpperReal,UpperReal);
    Real max(Real,UpperReal);
    Real max(UpperReal,Real);

    UpperReal min(UpperReal,UpperReal);

    UpperReal neg(LowerReal);
    LowerReal neg(UpperReal);
    UpperReal add(UpperReal,UpperReal);
    UpperReal sub(UpperReal,LowerReal);
    LowerReal sub(LowerReal,UpperReal);
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
    PositiveReal max(PositiveReal,PositiveReal);
    PositiveReal max(Real,PositiveReal);
    PositiveReal max(PositiveReal,Real);

    PositiveReal min(PositiveReal,PositiveReal);

    PositiveReal rec(PositiveReal);
    PositiveReal add(PositiveReal,PositiveReal);
    PositiveReal mul(PositiveReal,PositiveReal);
    PositiveReal div(PositiveReal,PositiveReal);
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
    PositiveLowerReal rec(PositiveUpperReal);
    PositiveUpperReal rec(PositiveLowerReal);
    PositiveLowerReal add(PositiveLowerReal,PositiveLowerReal);
    PositiveLowerReal mul(PositiveLowerReal,PositiveLowerReal);
    PositiveLowerReal div(PositiveLowerReal,PositiveUpperReal);
    PositiveUpperReal div(PositiveUpperReal,PositiveLowerReal);
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
    PositiveUpperReal rec(PositiveLowerReal);
    PositiveLowerReal rec(PositiveUpperReal);
    PositiveUpperReal add(PositiveUpperReal,PositiveUpperReal);
    PositiveUpperReal mul(PositiveUpperReal,PositiveUpperReal);
    PositiveUpperReal div(PositiveUpperReal,PositiveLowerReal);
    PositiveLowerReal div(PositiveLowerReal,PositiveUpperReal);
};

} // namespace Ariadne

#include "numeric/logical.h"

namespace Ariadne {

template<class N, EnableIf<IsIntegral<N>> =dummy> inline ValidatedKleenean operator==(const Real& x1, N n2) { return x1==Real(n2); }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline ValidatedKleenean operator!=(const Real& x1, N n2) { return x1!=Real(n2); }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline ValidatedKleenean operator<=(const Real& x1, N n2) { return x1<=Real(n2); }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline ValidatedKleenean operator>=(const Real& x1, N n2) { return x1>=Real(n2); }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline ValidatedKleenean operator< (const Real& x1, N n2) { return x1< Real(n2); }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline ValidatedKleenean operator> (const Real& x1, N n2) { return x1> Real(n2); }

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
