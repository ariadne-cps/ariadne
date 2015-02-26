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
template<> struct IsNumber<Real> : True { };

struct Accuracy { Nat _bits; Nat bits() const { return _bits; } TwoExp error() const; };

extern const Real pi;
extern const Real infinity;

//! \ingroup UserNumberSubModule
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
    typedef Effective Paradigm;
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

    explicit Real(ExactFloat64 x);

    operator Number<Effective>() const;

    // Extract floating-point properties
    UpperFloat64 upper() const;
    LowerFloat64 lower() const;
    ApproximateFloat64 approx() const;
    double get_d() const;

    // Extract arbitrarily accurate approximations
    BoundFloat64 operator() (Precision64 pr) const;
    BoundFloatMP operator() (PrecisionMP pr) const;
    BoundFloat64 get(Precision64 pr) const;
    BoundFloatMP get(PrecisionMP pr) const;
    BoundFloatMP evaluate(Accuracy acc) const;

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
    friend Real abs(Real);

    friend ErrorFloat64 mag(Real);

    friend NegSierpinski eq(Real,Real);
    friend Tribool lt(Real,Real);

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

    friend Number<Effective> operator+(Number<Effective>, Number<Effective>);
    friend Number<Effective> operator-(Number<Effective>, Number<Effective>);
    friend Number<Effective> operator*(Number<Effective>, Number<Effective>);
    friend Number<Effective> operator/(Number<Effective>, Number<Effective>);
  private:
    Real(std::int64_t n, Void*);
    Real(std::uint64_t m, Void*);
};

template<class M, EnableIf<And<IsIntegral<M>,IsUnsigned<M>>>> inline Real::Real(M m) : Real(std::uint64_t(m),nullptr) { };
template<class N, EnableIf<And<IsIntegral<N>,IsSigned<N>>>> inline Real::Real(N n) : Real(std::int64_t(n),nullptr) { };

Real pow(Real,Nat);
Real sqrt(Real);
Real exp(Real);
Real log(Real);
Real sin(Real);
Real cos(Real);
Real tan(Real);
Real atan(Real);



//! \ingroup UserNumberSubModule
//! \brief Computable lower real numbers defined by conversion to concrete floats.
class LowerReal
{
  private: public:
    SharedPointer<Real::Interface> _ptr;
  private: public:
    explicit LowerReal(SharedPointer<Real::Interface>);
  public:
    typedef EffectiveLower Paradigm;
    typedef LowerReal NumericType;
  public:
    LowerFloat64 operator() (Precision64 pr) const;
    LowerFloatMP operator() (PrecisionMP pr) const;
    LowerFloat64 get(Precision64 pr) const;
    LowerFloatMP get(PrecisionMP pr) const;
};

//! \ingroup UserNumberSubModule
//! \brief Computable lower real numbers defined by conversion to concrete floats.
class UpperReal
{
  private: public:
    SharedPointer<Real::Interface> _ptr;
  private: public:
    explicit UpperReal(SharedPointer<Real::Interface>);
  public:
    typedef EffectiveUpper Paradigm;
    typedef UpperReal NumericType;
  public:
    UpperFloat64 operator() (Precision64 pr) const;
    UpperFloatMP operator() (PrecisionMP pr) const;
    UpperFloat64 get(Precision64 pr) const;
    UpperFloatMP get(PrecisionMP pr) const;
};

} // namespace Ariadne

#include "numeric/logical.h"

namespace Ariadne {

template<class N, EnableIf<IsIntegral<N>> =dummy> inline Tribool operator==(const Real& x1, N n2) { return x1==Real(n2); }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline Tribool operator!=(const Real& x1, N n2) { return x1!=Real(n2); }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline Tribool operator<=(const Real& x1, N n2) { return x1<=Real(n2); }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline Tribool operator>=(const Real& x1, N n2) { return x1>=Real(n2); }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline Tribool operator< (const Real& x1, N n2) { return x1< Real(n2); }
template<class N, EnableIf<IsIntegral<N>> =dummy> inline Tribool operator> (const Real& x1, N n2) { return x1> Real(n2); }

/*
template<class D, EnableIf<IsFloatingPoint<D>> =dummy> inline auto
    operator==(Real r, D d) -> decltype(r==ApprxFloat64(d)) { return r==ApprxFloat64(d); }
template<class D, EnableIf<IsFloatingPoint<D>> =dummy> inline auto
    operator!=(Real r, D d) -> decltype(r!=ApprxFloat64(d)) { return r!=ApprxFloat64(d); }
template<class D, EnableIf<IsFloatingPoint<D>> =dummy> inline auto
    operator< (Real r, D d) -> decltype(r< ApprxFloat64(d)) { return r< ApprxFloat64(d); }
template<class D, EnableIf<IsFloatingPoint<D>> =dummy> inline auto
    operator> (Real r, D d) -> decltype(r> ApprxFloat64(d)) { return r> ApprxFloat64(d); }
template<class D, EnableIf<IsFloatingPoint<D>> =dummy> inline auto
    operator<=(Real r, D d) -> decltype(r<=ApprxFloat64(d)) { return r<=ApprxFloat64(d); }
template<class D, EnableIf<IsFloatingPoint<D>> =dummy> inline auto
    operator>=(Real r, D d) -> decltype(r>=ApprxFloat64(d)) { return r>=ApprxFloat64(d); }
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
