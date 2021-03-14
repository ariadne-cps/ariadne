/***************************************************************************
 *            numeric/reals.hpp
 *
 *  Copyright  2013-20  Pieter Collins
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

/*! \file numeric/reals.hpp
 *  \brief Real numbers, including lower, upper and naive types, and validated versions.
 */



#ifndef ARIADNE_REALS_HPP
#define ARIADNE_REALS_HPP


#include "utility/typedefs.hpp"
#include "utility/pointer.hpp"
#include "utility/handle.hpp"

#include "numeric/logical.decl.hpp"
#include "numeric/number.decl.hpp"
#include "numeric/float.decl.hpp"

#include "numeric/paradigm.hpp"
#include "numeric/arithmetic.hpp"
#include "numeric/sequence.hpp"

#include "numeric/dyadic.hpp"
#include "numeric/real.hpp"

namespace Ariadne {

//! A group
template<class G> class Group {
  public:
    G ineg() const; //!< \brief Inplace negation.
    friend G add(G,G); //!< \brief Addition
    friend G neg(G); //!< \brief Negation \relates Group
    friend OutputStream& operator<<(OutputStream& os, Group<G> const&); //!< \brief Write.
};

//! Symmetries
class Symmetry : public Group<Symmetry> {
  public:
    friend Symmetry sub(Symmetry,Symmetry); //!< \brief Subtraction
};

class Accuracy;

class Real;
class LowerReal;
class UpperReal;
class NaiveReal;
class PositiveReal;
class PositiveLowerReal;
class PositiveUpperReal;
class PositiveNaiveReal;
template<> struct IsNumber<Real> : True { };

class ValidatedReal;
class ValidatedUpperReal;
class ValidatedLowerReal;
class ApproximateReal;


class RealInterface;
using LowerRealInterface = RealInterface;
using UpperRealInterface = RealInterface;
using NaiveRealInterface = RealInterface;



//! \ingroup NumericModule
//! \brief Lower real number type \f$\R_<\f$.
//! \details An effectively computable <em>lower real</em> is a real number for which it is possible to compute arbitrarily accurate (dyadic) rational lower bounds;
//! equivalently, a bounded increasing sequence of (dyadic) rationals converging to the number.
//! However, the convergence rate is not known, and there may be arbitrarily large upward jumps in the sequence.
//!
//! Multiplication and division are unsupported, since given \f$x \geq \underline{x}\f$ and \f$y\geq\underline{y}\f$,
//! we cannot deduce \f$x \times y\geq\underline{x}\times\underline{y}\f$ if the bounds are negative.
//! However, multiplication and reciprocation of positive lower reals <em>is</em> supported.
//!
//! %Positive lower real numbers are a natural type for representing the measure of an open set,
//! or the separation between two points in a metric space.
//! \sa Real, UpperReal, NaiveReal
class LowerReal
    : public Handle<const LowerRealInterface>
    , public DirectedAbelian<LowerReal,UpperReal>
{
  public:
    typedef EffectiveTag Paradigm;
    typedef LowerReal NumericType;
  public:
    LowerReal(Real);
    explicit LowerReal(SharedPointer<const Interface>);
  public:
    //!@{
    //! \name Computation of rigorous validated lower bounds on the number
    ValidatedLowerReal compute(Effort eff) const;
    DyadicLowerBound compute_get(Effort eff) const;
    FloatDPLowerBound compute_get(Effort eff, DoublePrecision pr) const;
    FloatMPLowerBound compute_get(Effort eff, MultiplePrecision pr) const;
    //!@}
  public:
    //!@{
    //! \name Lattice operations
    friend PositiveNaiveReal abs(LowerReal const& r) = delete; //!< \em No absolute value operator!
    friend LowerReal min(LowerReal const& r1, LowerReal const& r2); //!< The mimimum of \a r1 and \a r2.
    friend LowerReal max(LowerReal const& r1, LowerReal const& r2); //!< The maximum of \a r1 and \a r2.
    //!@}

    //!@{
    //! \name Named arithmetical functions
    friend LowerReal pos(LowerReal const& r); //!< Identity \a +r.
    friend UpperReal neg(LowerReal const& r); //!< Negative \a -r.
    friend LowerReal neg(UpperReal const& r); //!< Negative \a -r.
    friend LowerReal hlf(LowerReal const& r); //!< Half \a r÷2.
    friend LowerReal add(LowerReal const& r1, LowerReal const& r2); //!< \brief Sum \a r1+r2.
    friend LowerReal sub(LowerReal const& r1, UpperReal const& r2); //!< \brief Difference \a r1-r2.
    friend UpperReal sub(UpperReal const& r1, LowerReal const& r2); //!< \brief Difference \a r1-r2.

    friend NaiveReal mul(LowerReal const& r1, LowerReal const& r2) = delete; //!< \brief \em No multiplication operator, since non-monotone!
    friend NaiveReal div(LowerReal const& r1, UpperReal const& r2) = delete; //!< \brief \em No division operator, since non-monotone!
    friend NaiveReal rec(LowerReal const& r) = delete; //!< \brief \em No reciprocal operation, since non-monotone!
    //!@}

    //!@{
    //! \name Named arithmetical functions on positive numbers
    friend LowerReal mul(LowerReal const& r1, PositiveReal const& r2); //!< Multiplication \a r1×r2 preserves monotonicity.
    friend LowerReal mul(PositiveReal const& r1, LowerReal const& r2); //!< Multiplication \a r1×r2 preserves monotonicity.
    friend LowerReal div(LowerReal const& r1, PositiveReal const& r2); //!< Division \a r1÷r2 preserves monotonicity.
    friend LowerReal div(PositiveReal const& r1, UpperReal const& r2); //!< Division \a r1÷r2 preserves monotonicity.
    friend UpperReal div(PositiveReal const& r1, LowerReal const& r2); //!< Division \a r1÷r2 preserves monotonicity.

    friend PositiveLowerReal add(PositiveLowerReal const& r1, PositiveLowerReal const& r2); //!< Addition \a r1+r2 preserves positivity.
    friend PositiveLowerReal mul(PositiveLowerReal const& r1, PositiveLowerReal const& r2); //!< Multiplication \a r1×r2 preserves positivity.
    friend PositiveLowerReal div(PositiveLowerReal const& r1, PositiveUpperReal const& r2); //!< Division \a r1÷r2 preserves positivity.
    friend PositiveUpperReal div(PositiveUpperReal const& r1, PositiveLowerReal const& r2); //!< Division \a r1÷r2 preserves positivity.
    friend PositiveLowerReal rec(PositiveUpperReal const& r); //!< Reciprocal \a 1/r preserves positivity.
    friend PositiveUpperReal rec(PositiveLowerReal const& r); //!< Reciprocal \a 1/r preserves positivity.
    //!@}

    //!@{
    //! \name Algebraic and transcendental functions
    friend PositiveLowerReal sqrt(PositiveLowerReal const& r); //!< The square root of \a r, √\a r. Requires \a r ≥ 0.
    friend PositiveLowerReal exp(LowerReal const& r); //!< The natural exponent of \a r, \em e<sup>r</sup>.
    friend LowerReal log(PositiveLowerReal const& r); //!< The natural logarithm of \a r. Requires \a r ≥ 0.
    friend LowerReal atan(LowerReal const& r); //!< The arc-tangent of \a r.
    friend PositiveLowerReal atan(PositiveLowerReal const& r); //!< Arc-tangent preserves positivity

    friend NaiveReal sin(LowerReal const& r); //!< No sine operation, since non-monotone!
    friend NaiveReal cos(LowerReal const& r); //!< No cosine operation, since non-monotone!
    friend NaiveReal tan(LowerReal const& r); //!< No tangent operation, since non-monotone!
    //!@}

    //!@{
    //! \name Limit operations.
    friend LowerReal limit(IncreasingSequence<Dyadic> const& qlbs);
        //!< Create a real number from an increasing sequence of dyadic lower bounds.
    //!@}

    //!@{
    //! \name Input/output operations
    friend OutputStream& operator<<(OutputStream& os, LowerReal const& r); //< Write to an output stream.
    //!@}

    //!@{
    //! \name Mixed operations
    friend Real min(LowerReal const& lr1, Real const& r2);
    friend Real min(Real const& r1, LowerReal const& lr2);
    //!@}
};

//! \ingroup NumericModule
//! \brief Upper real number type \f$\R_>\f$.
//! \details An effectively computable <em>upper real</em> is a real number for which it is possible to compute arbitrarily accurate (dyadic) rational upper bounds;
//! equivalently, a bounded decreasing sequence of (dyadic) rationals converging to the number.
//! However, the convergence rate is not known, and there may be arbitrarily large downward jumps in the sequence.
//!
//! Multiplication and division are unsupported, since given \f$x \leq \overline{x}\f$ and \f$y\leq\overline{y}\f$,
//! we cannot deduce \f$x \times y\leq\overline{x}\times\overline{y}\f$.
//! However, multiplication and reciprocation of positive upper reals <em>is</em> supported.
//!
//! %Positive upper real numbers are a natural type for representing upper bounds of a distance or norm,
//! including the radius of a ball, or the measure of an closed set.
//! \sa Real, LowerReal, NaiveReal
class UpperReal
    : public Handle<const UpperRealInterface>
    , public DirectedAbelian<UpperReal,LowerReal>
{
  public:
    typedef EffectiveTag Paradigm;
    typedef UpperReal NumericType;
  public:
    UpperReal(Real);
    explicit UpperReal(SharedPointer<const Interface>);
  public:
    //!@{
    //! \name Computation of rigorous validated upper bounds on the number
    ValidatedUpperReal compute(Effort eff) const;
    DyadicUpperBound compute_get(Effort eff) const;
    FloatDPUpperBound compute_get(Effort eff, DoublePrecision pr) const;
    FloatMPUpperBound compute_get(Effort eff, MultiplePrecision pr) const;
    //!@}

  public:
    //!@{
    //! \name Lattice operations
    friend PositiveNaiveReal abs(UpperReal const& r) = delete; //!< \em No absolute value operator!
    friend UpperReal min(UpperReal const& r1, UpperReal const& r2); //!< The mimimum of \a r1 and \a r2.
    friend UpperReal max(UpperReal const& r1, UpperReal const& r2); //!< The maximum of \a r1 and \a r2.
    //!@}

    //!@{
    //! \name Named arithmetical functions
    friend UpperReal pos(UpperReal const& r); //!< Identity \a +r.
    friend LowerReal neg(UpperReal const& r); //!< Negative \a -r.
    friend UpperReal neg(LowerReal const& r); //!< Negative \a -r.
    friend UpperReal hlf(UpperReal const& r); //!< Half \a r÷2.
    friend UpperReal add(UpperReal const& r1, UpperReal const& r2); //!< \brief Sum \a r1+r2.
    friend UpperReal sub(UpperReal const& r1, LowerReal const& r2); //!< \brief Difference \a r1-r2.
    friend LowerReal sub(LowerReal const& r1, UpperReal const& r2); //!< \brief Difference \a r1-r2.

    friend NaiveReal mul(UpperReal const& r1, UpperReal const& r2) = delete; //!< \brief \em No multiplication operator, since non-monotone!
    friend NaiveReal mul(UpperReal const& r1, LowerReal const& r2) = delete; //!< \brief \em No division operator, since non-monotone!
    friend NaiveReal rec(UpperReal const& r) = delete; //!< \brief \em No reciprocal operation, since non-monotone!
    //!@}

    //!@{
    //! \name Named arithmetical functions on positive numbers
    friend UpperReal mul(UpperReal const& r1, PositiveReal const& r2); //!< Multiplication \a r1×r2 preserves monotonicity.
    friend UpperReal mul(PositiveReal const& r1, UpperReal const& r2); //!< Multiplication \a r1×r2 preserves monotonicity.
    friend UpperReal div(UpperReal const& r1, PositiveReal const& r2); //!< Division \a r1÷r2 preserves monotonicity.
    friend UpperReal div(PositiveReal const& r1, LowerReal const& r2); //!< Division \a r1÷r2 preserves monotonicity.
    friend LowerReal div(PositiveReal const& r1, UpperReal const& r2); //!< Division \a r1÷r2 preserves monotonicity.

    friend PositiveUpperReal add(PositiveUpperReal const& r1, PositiveUpperReal const& r2); //!< Addition \a r1+r2 preserves positivity.
    friend PositiveUpperReal mul(PositiveUpperReal const& r1, PositiveUpperReal const& r2); //!< Multiplication \a r1×r2 preserves positivity.
    friend PositiveUpperReal div(PositiveUpperReal const& r1, PositiveLowerReal const& r2); //!< Division \a r1÷r2 preserves positivity.
    friend PositiveLowerReal div(PositiveLowerReal const& r1, PositiveUpperReal const& r2); //!< Division \a r1÷r2 preserves positivity.
    friend PositiveUpperReal rec(PositiveLowerReal const& r); //!< Reciprocal \a 1/r preserves positivity.
    friend PositiveLowerReal rec(PositiveUpperReal const& r); //!< Reciprocal \a 1/r preserves positivity.
    //!@}

    //!@{
    //! \name Algebraic and transcendental functions
    friend PositiveUpperReal sqrt(PositiveUpperReal const& r); //!< The square root of \a r, √\a r. Requires \a r ≥ 0.
    friend PositiveUpperReal exp(UpperReal const& r); //!< The natural exponent of \a r, \em e<sup>r</sup>.
    friend UpperReal log(PositiveUpperReal const& r); //!< The natural logarithm of \a r. Requires \a r ≥ 0.
    friend UpperReal atan(UpperReal const& r); //!< The arc-tangent of \a r.
    friend PositiveUpperReal atan(PositiveUpperReal const& r); //!< Arc-tangent preserves positivity
    friend NaiveReal sin(UpperReal const& r); //!< No sine operation, since non-monotone!
    friend NaiveReal cos(UpperReal const& r); //!< No cosine operation, since non-monotone!
    friend NaiveReal tan(UpperReal const& r); //!< No tangent operation, since non-monotone!
    //!@}

    //!@{
    //! \name Limit operations.
    friend UpperReal limit(DecreasingSequence<Dyadic> const& qubs);
        //!< Create a real number from an decreasing sequence of dyadic upper bounds.
    //!@}

    //!@{
    //! \name Input/output operations
    friend OutputStream& operator<<(OutputStream& os, UpperReal const& r); //< Write to an output stream.
    //!@}

    //!@{
    //! \name Mixed operations
    friend Real max(UpperReal const& ur1, Real const& r2);
    friend Real max(Real r1, UpperReal const& ur2);
    //!@}
};

//! \ingroup NumericModule
//! \brief %Real number type defined as limits of convergent sequences of (dyadic) rationals, without bounds on the convergence rate.
//! \details It is possible to compute arbitrarily accurate (dyadic) rational approximations, but the convergence rate is not known,
//! so at no point in the the computation can anything concrete be deduced about the value of the number.
//!
//! In principle useless for rigorous computation, but quickly-computed approximations to real numbers may be useful for preconditioning rigorous algorithms.
//! \sa Real, LowerReal, UpperReal
class NaiveReal
    : public Handle<const NaiveRealInterface>
    , public DeclareRealOperations<NaiveReal,PositiveNaiveReal>
    , public DeclareAnalyticFieldOperations<NaiveReal>
    , public DeclareLatticeOperations<NaiveReal,PositiveNaiveReal>
    , public DeclareComparisonOperations<NaiveReal,ApproximateKleenean>
    , public DefineFieldOperators<NaiveReal>
{
  public:
    typedef NaiveRealInterface Interface;
    typedef ApproximateTag Paradigm;
    typedef NaiveReal NumericType;
  public:
    //!@{
    //! \name Constructors
    NaiveReal(); //!< Default constructor yields the integer \c 0 as a real number.
    explicit NaiveReal(SharedPointer<const Interface>); //!< Construct from any class implementing the real number interface.
    //explicit NaiveReal(ConvergentSequence<DyadicBounds> const&); //!< Construct from a sequence of dyadic bounds converging to a singleton intersection.
    //!@}

    //!@{
    //! \name Conversion operators
    NaiveReal(Dbl d); //!< Construction from a builtin double-precision approximation.
    NaiveReal(Integer const& n); //!< Construct from an integer.
    NaiveReal(Dyadic const& d); //!< Construct from a dyadic number \a d=p/2<sup>q</sup> for integers \a p, \a q.
    NaiveReal(Decimal const& d); //!< Construct from a decimal number, given by its decimal expansion.
    NaiveReal(Rational const& q); //!< Construct from a rational number.

    NaiveReal(Real const& r); //!< Construct from an (effective) real number.
    NaiveReal(LowerReal const& lr); //!< Construct from a lower real number.
    NaiveReal(UpperReal const& ur); //!< Construct from an upper real number.

//    operator Number<Tag>() const; //!< Convert to an effective number for computations.
    //!@}

    //!@{
    //! \name Computation of approximations to the number.
    //! Compute an approximation with no guarantees on the error. The error converges to \a 0 as \a eff approaches infinity.
    //! The time taken should be roughly polynomial in \a eff.
    ApproximateReal compute(Effort eff) const;
    //!@}

    //!@{
    //! \name Standard arithmetic operators
    friend NaiveReal operator+(NaiveReal const& r); //!< Unary plus.
    friend NaiveReal operator-(NaiveReal const& r); //!< Unary minus.
    friend NaiveReal operator+(NaiveReal const& r1, NaiveReal const& r2); //!< Plus.
    friend NaiveReal operator-(NaiveReal const& r1, NaiveReal const& r2); //!< Minus.
    friend NaiveReal operator*(NaiveReal const& r1, NaiveReal const& r2); //!< Times.
    friend NaiveReal operator/(NaiveReal const& r1, NaiveReal const& r2); //!< Divides.
    friend NaiveReal& operator+=(NaiveReal& r1, NaiveReal const& r2); //!< Inplace plus.
    friend NaiveReal& operator-=(NaiveReal& r1, NaiveReal const& r2); //!< Inplace minus.
    friend NaiveReal& operator*=(NaiveReal& r1, NaiveReal const& r2); //!< Inplace times.
    friend NaiveReal& operator/=(NaiveReal& r1, NaiveReal const& r2); //!< Inplace divides.
    //!@}

    //!@{
    //! \name Named arithmetical functions
    friend NaiveReal pos(NaiveReal const& r); //!< Identity \a +r.
    friend NaiveReal neg(NaiveReal const& r); //!< Negative \a -r.
    friend NaiveReal hlf(NaiveReal const& r); //!< Half \a r÷2.
    friend NaiveReal sqr(NaiveReal const& r); //!< Square \a r<sup>2</sup>.
    friend NaiveReal rec(NaiveReal const& r); //!< Reciprocal \a 1/r.
    friend NaiveReal add(NaiveReal const& r1, NaiveReal const& r2); //!< \brief Sum \a r1+r2.
    friend NaiveReal sub(NaiveReal const& r1, NaiveReal const& r2); //!< \brief Difference \a r1-r2.
    friend NaiveReal mul(NaiveReal const& r1, NaiveReal const& r2); //!< \brief Product \a r1×r2.
    friend NaiveReal div(NaiveReal const& r1, NaiveReal const& r2); //!< \brief Quotient \a r1÷r2.
    friend NaiveReal fma(NaiveReal const& r1, NaiveReal const& r2, NaiveReal const& r3); //!< \brief Fused multiply-and-add \a r1×r2+r3.
    friend NaiveReal pow(NaiveReal const& r, Int n); //!< \brief Power \a r<sup>n</sup>.
    //!@}

    //!@{
    //! \name Algebraic and transcendental functions
    friend NaiveReal sqrt(NaiveReal const& r); //!< The square root of \a r, √\a r. Requires \a r ≥ 0.
    friend NaiveReal exp(NaiveReal const& r); //!< The natural exponent of \a r, \em e<sup>r</sup>.
    friend NaiveReal log(NaiveReal const& r); //!< The natural logarithm of \a r. Requires \a r ≥ 0.
    friend NaiveReal sin(NaiveReal const& r); //!< The sine of \a r.
    friend NaiveReal cos(NaiveReal const& r); //!< The cosine of \a r.
    friend NaiveReal tan(NaiveReal const& r); //!< The tangent of \a r, sin(\a r)/cos(\a r).
    friend NaiveReal atan(NaiveReal const& r); //!< The arc-tangent of \a r.
    //!@}

    //!@{
    //! \name Lattice operations
    friend PositiveNaiveReal abs(NaiveReal const&); //!< Absolute value \a |r|.
    friend NaiveReal min(NaiveReal const& r1, NaiveReal const& r2); //!< The mimimum of \a r1 and \a r2.
    friend NaiveReal max(NaiveReal const& r1, NaiveReal const& r2); //!< The maximum of \a r1 and \a r2.
    //!@}

    //!@{
    //! Operations based on the metric structure.
    friend PositiveNaiveReal dist(NaiveReal const& r1, NaiveReal const& r2);
        //< The distance |\a r <sub>1</sub>-\a r <sub>2</sub>| between \a r<sub>1</sub> and \a r<sub>2</sub>.
    //!@}

    //!@{
    //! \name Comparison operations and operators.
    friend ApproximateKleenean leq(NaiveReal const& r1, NaiveReal const& r2);
        //!< Returns \c likely if \a r1<r2, \c unlikely if \a r1\>r2, and \c indetermiate if \a r1==r2.
    //!@}

    //!@{
    //! \name Limit operations.
    friend NaiveReal limit(ConvergentSequence<Dyadic> const& qas);
        //!< Create a real number from an converging sequence of dyadic approximants.
    //!@}

    //!@{
    //! \name Input/output operations
    friend OutputStream& operator<<(OutputStream& os, NaiveReal const& r); //< Write to an output stream.
    //!@}
};

//! \ingroup UserNumericTypeSubModule
//! \brief Computable lower real numbers defined by conversion to concrete floats.
class PositiveLowerReal : public LowerReal, public DirectedSemiRing<PositiveLowerReal,PositiveUpperReal>
{
  public:
    PositiveLowerReal(PositiveReal r) : LowerReal(r) { }
    explicit PositiveLowerReal(LowerReal r) : LowerReal(r) { }
    PositiveLowerBound<Dyadic> compute_get(Effort) const;
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
class PositiveUpperReal : public UpperReal, public DirectedSemiRing<PositiveUpperReal,PositiveLowerReal>
{
  public:
    PositiveUpperReal(PositiveReal r) : UpperReal(r) { }
    explicit PositiveUpperReal(UpperReal r) : UpperReal(r) { }
    PositiveUpperBound<Dyadic> compute_get(Effort) const;
  public:
    PositiveUpperReal rec(PositiveLowerReal const&);
    PositiveLowerReal rec(PositiveUpperReal const&);
    PositiveUpperReal add(PositiveUpperReal const&, PositiveUpperReal const&);
    PositiveUpperReal mul(PositiveUpperReal const&, PositiveUpperReal const&);
    PositiveUpperReal div(PositiveUpperReal const&, PositiveLowerReal const&);
    PositiveLowerReal div(PositiveLowerReal const&, PositiveUpperReal const&);
};

PositiveUpperReal rec(PositiveLowerReal plr);
PositiveLowerReal rec(PositiveUpperReal pur);
PositiveLowerReal add(PositiveLowerReal plr1, PositiveLowerReal plr2);
PositiveUpperReal add(PositiveUpperReal pur1, PositiveUpperReal pur2);
PositiveLowerReal mul(PositiveLowerReal plr1, PositiveLowerReal plr2);
PositiveUpperReal mul(PositiveUpperReal pur1, PositiveUpperReal pur2);
PositiveLowerReal div(PositiveLowerReal plr1, PositiveUpperReal pur2);
PositiveUpperReal div(PositiveUpperReal pur1, PositiveLowerReal plr2);



class ValidatedRealInterface;

//! \ingroup NumericModule
//! \brief A generic class representing rigorous bounds on a real number.
//! \see Real
class ValidatedReal
    : public Handle<ValidatedRealInterface>
{
  public:
    typedef ValidatedRealInterface Interface;

    ValidatedReal(DyadicBounds const&);
    operator DyadicBounds() const;

    Dyadic value();
    Dyadic error();
    Dyadic lower();
    Dyadic upper();

    //! \brief Get dyadic bounds for the number.
    //! \details It is not always clear how this function can be implemented,
    //!   and it may be removed or modified in future versions.
    DyadicBounds get() const;
    //! \brief Get the bounds for the number, representing in double precision.
    FloatDPBounds get(DoublePrecision) const;
    //! \brief Get the bounds for the number, representing to precision \a pr.
    //! Note that increasing the accuracy typically does not yield arbitrarily
    //! tight bounds, as the object is already an approximation to finite accuracy.
    FloatMPBounds get(MultiplePrecision pr) const;
    //! \brief Write to an output stream.
    friend OutputStream& operator<<(OutputStream&, ValidatedReal const&);

    friend ValidatedKleenean sgn(DyadicBounds const& r);
    friend ValidatedKleenean operator==(DyadicBounds const& vr1, DyadicBounds const& vr2);
    friend ValidatedKleenean operator!=(DyadicBounds const& vr1, DyadicBounds const& vr2);
    friend ValidatedKleenean operator<=(DyadicBounds const& vr1, DyadicBounds const& vr2);
    friend ValidatedKleenean operator>=(DyadicBounds const& vr1, DyadicBounds const& vr2);
    friend ValidatedKleenean operator< (DyadicBounds const& vr1, DyadicBounds const& vr2);
    friend ValidatedKleenean operator> (DyadicBounds const& vr1, DyadicBounds const& vr2);
};

class ValidatedLowerReal
    : public Handle<ValidatedRealInterface>
{
  public:
    ValidatedLowerReal(DyadicLowerBound const&);
    ValidatedLowerReal(ValidatedReal const&);
    DyadicLowerBound get() const;
    FloatDPLowerBound get(DoublePrecision) const;
    FloatMPLowerBound get(MultiplePrecision) const;
    friend OutputStream& operator<<(OutputStream&, ValidatedLowerReal const&);
};

class ValidatedUpperReal
    : public Handle<ValidatedRealInterface>
{
  public:
    ValidatedUpperReal(DyadicUpperBound const&);
    ValidatedUpperReal(ValidatedReal const&);
    DyadicUpperBound get() const;
    FloatDPUpperBound get(DoublePrecision) const;
    FloatMPUpperBound get(MultiplePrecision) const;
    friend OutputStream& operator<<(OutputStream&, ValidatedUpperReal const&);
};

class ApproximateRealInterface;

//! \ingroup NumericModule
//! \brief A generic class representing an approximation to a real number
//! of unknown accuracy.
//! \see Real, ValidatedReal
class ApproximateReal
    : public Handle<ApproximateRealInterface>
{
  public:
    typedef ApproximateRealInterface Interface;

    ApproximateReal(DyadicApproximation const&);
    //! \brief Get a dyadic approximation to the number.
    DyadicApproximation get() const;
    //! \brief Get the approximation to the number, converting to double-precision.
    FloatDPApproximation get(DoublePrecision) const;
    //! \brief Get the approximation to the number, converting to a representation
    //! with precision \a pr.
    //! Note that increasing the precision typically does not yield arbitrarily
    //! accurate approximation, as the object already has some error.
    FloatMPApproximation get(MultiplePrecision pr) const;
    //! \brief Write to an output stream.
    friend OutputStream& operator<<(OutputStream&, ApproximateReal const&);
};


} // namespace Ariadne

#endif
