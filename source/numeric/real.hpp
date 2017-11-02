/***************************************************************************
 *            numeric/real.hpp
 *
 *  Copyright 2013-17  Pieter Collins
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

/*! \file numeric/real.hpp
 *  \brief
 */



#ifndef ARIADNE_REAL_HPP
#define ARIADNE_REAL_HPP

#include "utility/typedefs.hpp"
#include "utility/pointer.hpp"

#include "numeric/logical.decl.hpp"
#include "numeric/number.decl.hpp"
#include "numeric/float.decl.hpp"

#include "paradigm.hpp"
#include "numeric/arithmetic.hpp"
#include "numeric/sequence.hpp"

namespace Ariadne {

class Real;
class LowerReal;
class UpperReal;
class PositiveReal;
class PositiveLowerReal;
class PositiveUpperReal;
template<> struct IsNumericType<Real> : True { };

class ValidatedReal;

//! \ingroup NumericModule
//! \brief The accuracy of a computation of a value in a metric space. The result must be correct to within 2<sup>-acc</sup> bits.
struct Accuracy {
    Nat _bits;
    Accuracy(Nat bits) : _bits(bits) { }
    Nat bits() const { return _bits; }
    TwoExp error() const;
    friend OutputStream& operator<<(OutputStream& os, Accuracy acc) { return os << "Accuracy("<<acc.bits()<<")"; }
};

extern const Real pi;
extern const Real infinity;

class RealInterface;

//! \ingroup NumericModule
//! \brief %Real numbers supporting elementary functions, comparisons and limits.
class Real
    : public DeclareRealOperations<Real,PositiveReal>
    , public DeclareAnalyticFieldOperations<Real>
    , public DeclareLatticeOperations<Real,PositiveReal>
    , public DeclareComparisonOperations<Real,Kleenean,NegatedSierpinskian>
    , public DefineFieldOperators<Real>
{
  private: public:
    using Interface = RealInterface;
    using InterfaceType = RealInterface;
  private: public:
    SharedPointer<Interface> _ptr;
  private:
    explicit Real(double,double,double);
  public:
    typedef EffectiveTag Paradigm;
    typedef Real NumericType;
  public:
    //@{
    //! \name Constructors
    Real(); //!< Default constructor yields the integer \c 0 as a real number.
    explicit Real(SharedPointer<Real::InterfaceType>); //!< Construct from any class implementing the real number interface.
    explicit Real(ConvergentSequence<DyadicBounds> const&); //!< Construct from a sequence of dyadic bounds converging to a singleton intersection.
    explicit Real(StrongCauchySequence<Dyadic> const&); //!< Construct from a strongly convergent sequence of dyadic numbers.
    //@}

    //@{
    //! \name Conversion operators
    explicit Real(double); //!< Unsafe construction from a double is DEPRECATED

#ifdef DOXYGEN
    template<class N, EnableIf<IsBuiltinIntegral<N>> = dummy> Real(N n); //!< Convert from a builtin integer (but not a float).
#else
    template<class M, EnableIf<And<IsBuiltinIntegral<M>,IsBuiltinUnsigned<M>>> = dummy> Real(M m);
    template<class N, EnableIf<And<IsBuiltinIntegral<N>,IsBuiltinSigned<N>>> = dummy> Real(N n);
#endif
    Real(ExactDouble d); //!< Construct from a double-precision value representing a number exactly.
    Real(Integer const& n); //!< Construct from an integer.
    Real(Dyadic const& d); //!< Construct from a dyadic number \a d=p/2<sup>q</sup> for integers \a p, \a q.
    Real(Decimal const& d); //!< Construct from a decimal number, given by its decimal expansion.
    Real(Rational const& q); //!< Construct from a rational number.

    explicit Real(FloatDPValue x); //!< DEPRECATED
    explicit Real(EffectiveNumber r); //!< DEPRECATED

    operator Number<EffectiveTag>() const; //!< Convert to an effective number for computations.
    //@}

    //@{
    //! \name Convert to real number types describing less information.
    UpperReal upper() const; //!< A real number allowing only computation of upper bounds.
    LowerReal lower() const; //!< A real number allowing only computation of lower bounds.
    //FloatDPApproximation approx() const; //!< A real number allowing only computation of approximations with no guarantees on the accuracy.
    //@}

    //@{
    //! \name Computation of rigorous validated bounds on the number.

    //! Compute a concrete approximation with an error of at most \a 2<sup>-acc</sup> i.e. \a acc binary digits.
    ValidatedReal compute(Accuracy acc) const;
    //! Compute a concrete approximation with a bound on the error. The error bound converges to \a 0 as \a eff approaches infinity.
    //! The time taken should be roughly polynomial in \a eff.
    ValidatedReal compute(Effort eff) const;
    //! Compute a concrete approximation using double-precision.
    FloatDPBounds get(DoublePrecision pr) const;
    //! Compute a concrete approximation using the given precision.
    FloatMPBounds get(MultiplePrecision pr) const;
    //@}

    //@{
    //! \name Named arithmetical functions
    friend Real pos(Real const& r); //!< Identity \a +r.
    friend Real neg(Real const& r); //!< Negative \a -r.
    friend Real hlf(Real const& r); //!< Half \a r÷2.
    friend Real sqr(Real const& r); //!< Square \a r<sup>2</sup>.
    friend Real rec(Real const& r); //!< Reciprocal \a 1/r.
    friend Real add(Real const& r1, Real const& r2); //!< \brief Sum \a r1+r2.
    friend Real sub(Real const& r1, Real const& r2); //!< \brief Difference \a r1-r2.
    friend Real mul(Real const& r1, Real const& r2); //!< \brief Product \a r1×r2.
    friend Real div(Real const& r1, Real const& r2); //!< \brief Quotient \a r1÷r2.
    friend Real fma(Real const& r1, Real const& r2, Real const& r3); //!< \brief Fused multiply-and-add \a r1×r2+r3.
    friend Real pow(Real const& r, Int n); //!< \brief Power \a r<sup>n</sup>.
    //@}

    //@{
    //! \name Algebraic and transcendental functions
    friend Real sqrt(Real const& r); //!< The square root of \a r, √\a r. Requires \a r ≥ 0.
    friend Real exp(Real const& r); //!< The natural exponent of \a r, \em e<sup>r</sup>.
    friend Real log(Real const& r); //!< The natural logarithm of \a r. Requires \a r ≥ 0.
    friend Real atan(Real const& r); //!< The arc-tangent of \a r.
    friend Real sin(Real const& r); //!< The sine of \a r.
    friend Real cos(Real const& r); //!< The cosine of \a r.
    friend Real tan(Real const& r); //!< The tangent of \a r, sin(\a r)/cos(\a r) \f$.
    friend Real stan(Real const& r); //!< The arc-tangent of \a r.
    //@}

    //@{
    //! \name Standard arithmetic operators
    friend Real operator+(Real const& r); //!< Unary plus.
    friend Real operator-(Real const& r); //!< Unary minus.
    friend Real operator+(Real const& r1, Real const& r2); //!< Plus.
    friend Real operator-(Real const& r1, Real const& r2); //!< Minus.
    friend Real operator*(Real const& r1, Real const& r2); //!< Times.
    friend Real operator/(Real const& r1, Real const& r2); //!< Divides.
    friend Real& operator+=(Real& r1, Real const& r2); //!< Inplace plus.
    friend Real& operator-=(Real& r1, Real const& r2); //!< Inplace minus.
    friend Real& operator*=(Real& r1, Real const& r2); //!< Inplace times.
    friend Real& operator/=(Real& r1, Real const& r2); //!< Inplace divides.
    //@}

    //@{
    //! \name Lattice operations
    friend PositiveReal abs(Real const&); //!< Absolute value \a |r|.
    friend Real min(Real const& r1, Real const& r2); //!< The mimimum of \a r1 and \a r2.
    friend Real max(Real const& r1, Real const& r2); //!< The maximum of \a r1 and \a r2.
    //@}

    //@{
    //! Operations based on the metric structure.
    friend PositiveReal dist(Real const& r1, Real const& r2);
        //< The distance |\a r <sub>1</sub>-\a r <sub>2</sub>| between \a r<sub>1</sub> and \a r<sub>2</sub>.
    friend PositiveUpperReal mag(Real const& r);
        //!< An over-approximation to the absolute value of \a r.
    friend FloatDPError mag(Real const&, DoublePrecision);
    //@}

    //@{
    //! \name Comparison operations and operators.
    friend Kleenean leq(Real const& r1, Real const& r2); //!< Returns \c true if \a r1<r2, \c false if \a r1\>r2, and \c indetermiate if \a r1==r2.
    friend NegatedSierpinskian operator==(Real const& r1, Real const& r2); //!< Equality is undecidable and may only robustly be falsified.
    friend Sierpinskian operator!=(Real const& r1, Real const& r2); //!< Inequality is undecidable and may only robustly be verified.
    friend Kleenean operator<=(Real const& r1, Real const& r2); //!< Comparison \c leq.
    friend Kleenean operator>=(Real const& r1, Real const& r2); //!< Comparison .
    friend Boolean operator>(Real const& r1, Pair<Real,Real> lu); //!< Given \a l<u, returns \a true if r>l and \a false if r<u.
    //@}

    //@{
    //! \name Limit operations.
    friend Real limit(ConvergentSequence<DyadicBounds> const& qbs);
        //!< Create a real number from a sequence of dyadic bounds whose width converges to 0.
    friend Real limit(StrongCauchySequence<Dyadic> const& qs);
        //!< Create a real number from a sequence of dyadic numbers \em q<sub>n</sub> for which
        //!< <em>|q<sub>n<sub>1</sub></sub>-q<sub>n<sub>2</sub></sub>|≤<em>2<sup>-min</sup><sup>(n<sub>1</sub>,n<sub>2</sub>)</sup></em>
    friend Real limit(StrongCauchySequence<Real> const& rs);
        //!< The limit of a sequence of real numbers \em r<sub>n</sub> for which
        //!< <em>|r<sub>n<sub>1</sub></sub>-r<sub>n<sub>2</sub></sub>|≤<em>2<sup>-min</sup><sup>(n<sub>1</sub>,n<sub>2</sub>)</sup></em>
    //@}


    //@{
    //! \name Operations on the representation.
    friend Bool same(Real const&, Real const&); //!< Test equivalence of representation.
    //@}

    //@{
    //! \name Input/output operations
    friend OutputStream& operator<<(OutputStream& os, Real const& r); //< Write to an output stream.
    //@}
    double get_d() const;
/*
    friend Number<EffectiveTag> operator+(Number<EffectiveTag>, Number<EffectiveTag>);
    friend Number<EffectiveTag> operator-(Number<EffectiveTag>, Number<EffectiveTag>);
    friend Number<EffectiveTag> operator*(Number<EffectiveTag>, Number<EffectiveTag>);
    friend Number<EffectiveTag> operator/(Number<EffectiveTag>, Number<EffectiveTag>);

    template<class N, EnableIf<IsBuiltinIntegral<N>> =dummy> friend inline decltype(auto) operator==(const Real& x1, N n2) { return x1==Real(n2); }
    template<class N, EnableIf<IsBuiltinIntegral<N>> =dummy> friend inline decltype(auto) operator!=(const Real& x1, N n2) { return x1!=Real(n2); }
    template<class N, EnableIf<IsBuiltinIntegral<N>> =dummy> friend inline decltype(auto) operator<=(const Real& x1, N n2) { return x1<=Real(n2); }
    template<class N, EnableIf<IsBuiltinIntegral<N>> =dummy> friend inline decltype(auto) operator>=(const Real& x1, N n2) { return x1>=Real(n2); }
    template<class N, EnableIf<IsBuiltinIntegral<N>> =dummy> friend inline decltype(auto) operator< (const Real& x1, N n2) { return x1< Real(n2); }
    template<class N, EnableIf<IsBuiltinIntegral<N>> =dummy> friend inline decltype(auto) operator> (const Real& x1, N n2) { return x1> Real(n2); }
*/

  private:
    Real(std::int64_t n, Void*);
    Real(std::uint64_t m, Void*);
};

template<class M, EnableIf<And<IsBuiltinIntegral<M>,IsBuiltinUnsigned<M>>>> inline Real::Real(M m) : Real(std::uint64_t(m),nullptr) { };
template<class N, EnableIf<And<IsBuiltinIntegral<N>,IsBuiltinSigned<N>>>> inline Real::Real(N n) : Real(std::int64_t(n),nullptr) { };


//! \ingroup NumericModule
//! \brief Computable lower real numbers defined by conversion to concrete floats.
class LowerReal
    : public DirectedRing<LowerReal,UpperReal,PositiveLowerReal>
    , public DirectedRing<UpperReal,LowerReal,PositiveUpperReal>
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
    FloatDPLowerBound operator() (DoublePrecision pr) const;
    FloatMPLowerBound operator() (MultiplePrecision pr) const;
    FloatDPLowerBound get(DoublePrecision pr) const;
    FloatMPLowerBound get(MultiplePrecision pr) const;
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

//! \ingroup NumericModule
//! \brief Computable lower real numbers defined by conversion to concrete floats.
class UpperReal
    : public DirectedRing<UpperReal,LowerReal,PositiveUpperReal>
    , public DirectedRing<LowerReal,UpperReal,PositiveLowerReal>
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
    FloatDPUpperBound operator() (DoublePrecision pr) const;
    FloatMPUpperBound operator() (MultiplePrecision pr) const;
    FloatDPUpperBound get(DoublePrecision pr) const;
    FloatMPUpperBound get(MultiplePrecision pr) const;
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
    PositiveFloatDPBounds get(DoublePrecision pr) const;
    PositiveFloatMPBounds get(MultiplePrecision pr) const;
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
    PositiveFloatDPLowerBound get(DoublePrecision pr) const;
    PositiveFloatMPLowerBound get(MultiplePrecision pr) const;
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
    PositiveFloatDPUpperBound get(DoublePrecision pr) const;
    PositiveFloatMPUpperBound get(MultiplePrecision pr) const;
  public:
    PositiveUpperReal rec(PositiveLowerReal const&);
    PositiveLowerReal rec(PositiveUpperReal const&);
    PositiveUpperReal add(PositiveUpperReal const&, PositiveUpperReal const&);
    PositiveUpperReal mul(PositiveUpperReal const&, PositiveUpperReal const&);
    PositiveUpperReal div(PositiveUpperReal const&, PositiveLowerReal const&);
    PositiveLowerReal div(PositiveLowerReal const&, PositiveUpperReal const&);
};

PositiveReal cast_positive(Real const& x);


class ValidatedRealInterface;

class ValidatedReal {
    SharedPointer<ValidatedRealInterface> _ptr;
  public:
    ValidatedReal(DyadicBounds const&);
    Dyadic value();
    Dyadic error();
    Dyadic lower();
    Dyadic upper();
    DyadicBounds get() const;
    FloatDPBounds get(DoublePrecision) const;
    FloatMPBounds get(MultiplePrecision) const;
    friend OutputStream& operator<<(OutputStream&, ValidatedReal const&);
};

} // namespace Ariadne

#include "numeric/logical.hpp"

namespace Ariadne {

/*
template<class D, EnableIf<IsBuiltinFloatingPoint<D>> =dummy> inline auto
    operator==(Real r, D d) -> decltype(r==FloatDPApproximation(d)) { return r==FloatDPApproximation(d); }
template<class D, EnableIf<IsBuiltinFloatingPoint<D>> =dummy> inline auto
    operator!=(Real r, D d) -> decltype(r!=FloatDPApproximation(d)) { return r!=FloatDPApproximation(d); }
template<class D, EnableIf<IsBuiltinFloatingPoint<D>> =dummy> inline auto
    operator< (Real r, D d) -> decltype(r< FloatDPApproximation(d)) { return r< FloatDPApproximation(d); }
template<class D, EnableIf<IsBuiltinFloatingPoint<D>> =dummy> inline auto
    operator> (Real r, D d) -> decltype(r> FloatDPApproximation(d)) { return r> FloatDPApproximation(d); }
template<class D, EnableIf<IsBuiltinFloatingPoint<D>> =dummy> inline auto
    operator<=(Real r, D d) -> decltype(r<=FloatDPApproximation(d)) { return r<=FloatDPApproximation(d); }
template<class D, EnableIf<IsBuiltinFloatingPoint<D>> =dummy> inline auto
    operator>=(Real r, D d) -> decltype(r>=FloatDPApproximation(d)) { return r>=FloatDPApproximation(d); }
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
