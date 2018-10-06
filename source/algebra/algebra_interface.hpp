/***************************************************************************
 *            algebra_interface.hpp
 *
 *  Copyright 2010-17  Pieter Collins
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

/*! \file algebra_interface.hpp
 *  \brief Interface for function algebras.
 */

#ifndef ARIADNE_ALGEBRA_INTERFACE_HPP
#define ARIADNE_ALGEBRA_INTERFACE_HPP

#include <iosfwd>
#include <iostream>

#include "../utility/writable.hpp"
#include "../numeric/numeric.hpp"
#include "../utility/pointer.hpp"
#include "../numeric/operators.hpp"

namespace Ariadne {

template<class X> class AlgebraInterface;
template<class X> class Algebra;
template<class X> class NormedAlgebra;
template<class X> class GradedAlgebra;
template<class X> class SymbolicAlgebra;

template<class X> struct AlgebraTraits;

template<> struct AlgebraTraits<ApproximateNumber> {
    typedef ApproximateNumber ValueType;
    typedef Interval<ApproximateNumber> RangeType;
    typedef Positive<ApproximateNumber> NormType;
    typedef ApproximateNumber NumericType;
};

template<> struct AlgebraTraits<ValidatedNumber> {
    typedef ExactNumber ValueType;
    typedef Interval<ValidatedUpperNumber> RangeType;
    typedef Positive<ValidatedNumber> NormType;
    typedef ValidatedNumber NumericType;
};

template<> struct AlgebraTraits<EffectiveNumber> {
    typedef EffectiveNumber ValueType;
    typedef Interval<EffectiveUpperNumber> RangeType;
    typedef Positive<EffectiveNumber> NormType;
    typedef EffectiveNumber NumericType;
};

template<class F> struct AlgebraTraits<Approximation<F>> {
    typedef Approximation<F> ValueType;
    typedef Interval<Approximation<F>> RangeType;
    typedef PositiveApproximation<F> NormType;
    typedef Approximation<F> NumericType;
};

template<class F> struct AlgebraTraits<Bounds<F>> {
    typedef Value<F> ValueType;
    typedef Interval<UpperBound<F>> RangeType;
    typedef FloatDPError NormType;
    typedef Bounds<F> NumericType;
};

template<> struct AlgebraTraits<Real> {
    typedef FloatDPValue ValueType;
    typedef Interval<FloatDPUpperBound> RangeType;
    typedef FloatDPError NormType;
    typedef Real NumericType;
};


//! \brief Interface for a unital algebra over a field \a X.
template<class X> class AlgebraInterface
    : public virtual WritableInterface
{
    friend class Algebra<X>;
  public:
    typedef X NumericType;
  public:
    //! \brief Virtual destructor.
    virtual ~AlgebraInterface<X>() = default;
    //! \brief Create a dynamically-allocated copy.
    virtual AlgebraInterface<X>* _create_copy() const = 0;
    //! \brief Create the zero element in the same algebra as the current object.
    virtual AlgebraInterface<X>* _create_zero() const = 0;
    //! \brief Create a constant multiple of the unit element in the same algebra as the current object.
    virtual AlgebraInterface<X>* _create_constant(X const& c) const = 0;

    //! \brief Add a constant numerical scalar \c r+=c .
    virtual Void _iadd(const X& c) = 0;
    //! \brief Multiply by a numerical scalar \c r*=c .
    virtual Void _imul(const X& c) = 0;
    //! \brief Fused scalar multiply and add \c r+=x1*x2 .
    virtual Void _isma(const X& c1, const AlgebraInterface<X>& x2) = 0;
    //! \brief Fused multiply and add \c r+=x1*x2 .
    virtual Void _ifma(const AlgebraInterface<X>& x1, const AlgebraInterface<X>& x2) = 0;

    //! \brief Negation (unary minus).
    virtual AlgebraInterface<X>* _neg() const = 0;
    //! \brief Add another algebra element.
    virtual AlgebraInterface<X>* _add(const AlgebraInterface<X>& x) const = 0;
    //! \brief Subract another algebra element.
    virtual AlgebraInterface<X>* _sub(const AlgebraInterface<X>& x) const = 0;
    //! \brief Add multiply with another algebra element \c r+=c*x .
    virtual AlgebraInterface<X>* _mul(const AlgebraInterface<X>& x) const = 0;

    virtual AlgebraInterface<X>* _add(const X& c) const = 0;
    virtual AlgebraInterface<X>* _sub(const X& c) const = 0;
    virtual AlgebraInterface<X>* _mul(const X& c) const = 0;
    virtual AlgebraInterface<X>* _div(const X& c) const = 0;
    virtual AlgebraInterface<X>* _radd(const X& c) const = 0;
    virtual AlgebraInterface<X>* _rsub(const X& c) const = 0;
    virtual AlgebraInterface<X>* _rmul(const X& c) const = 0;

    virtual AlgebraInterface<X>* _pow(Nat m) const = 0;

    virtual AlgebraInterface<X>* _apply(Neg) const = 0;
    virtual AlgebraInterface<X>* _apply(Add,AlgebraInterface<X>const&) const = 0;
    virtual AlgebraInterface<X>* _apply(Sub,AlgebraInterface<X>const&) const = 0;
    virtual AlgebraInterface<X>* _apply(Mul,AlgebraInterface<X>const&) const = 0;
    virtual AlgebraInterface<X>* _apply(Add,X const&) const = 0;
    virtual AlgebraInterface<X>* _apply(Mul,X const&) const = 0;
};

//! \brief Interface for a normed unital algebra over a field \a X.
template<class X> class TranscendentalAlgebraInterface
    : public virtual AlgebraInterface<X>
{
  public:
    typedef typename AlgebraTraits<X>::ValueType ValueType;
    typedef typename AlgebraTraits<X>::NormType NormType;
    typedef typename AlgebraTraits<X>::RangeType RangeType;
  public:
    // Overrides for AlgebraInterface operations
    virtual TranscendentalAlgebraInterface<X>* _create_copy() const = 0;
    virtual TranscendentalAlgebraInterface<X>* _create_zero() const = 0;
    virtual TranscendentalAlgebraInterface<X>* _create_constant(X const& c) const = 0;

    using AlgebraInterface<X>::_apply;
    virtual TranscendentalAlgebraInterface<X>* _apply(Rec) const;
    virtual TranscendentalAlgebraInterface<X>* _apply(Sqrt) const;
    virtual TranscendentalAlgebraInterface<X>* _apply(Exp) const;
    virtual TranscendentalAlgebraInterface<X>* _apply(Log) const;
    virtual TranscendentalAlgebraInterface<X>* _apply(Sin) const;
    virtual TranscendentalAlgebraInterface<X>* _apply(Cos) const;
    virtual TranscendentalAlgebraInterface<X>* _apply(Tan) const;
    virtual TranscendentalAlgebraInterface<X>* _apply(Atan) const;
};

//! \brief Interface for a normed unital algebra over a field \a X.
template<class X> class NormedAlgebraInterface
    : public virtual AlgebraInterface<X>
{
  public:
    typedef typename AlgebraTraits<X>::ValueType ValueType;
    typedef typename AlgebraTraits<X>::NormType NormType;
    typedef typename AlgebraTraits<X>::RangeType RangeType;
  public:
    NormedAlgebra<X> create() const;
    NormedAlgebra<X> clone() const;
  public:
    // Overrides for AlgebraInterface operations
    virtual NormedAlgebraInterface<X>* _create_copy() const = 0;
    virtual NormedAlgebraInterface<X>* _create_zero() const = 0;
    virtual NormedAlgebraInterface<X>* _create_constant(X const& c) const = 0;
    virtual NormedAlgebraInterface<X>* _create_ball(ErrorType const& r) const = 0;

    //! \brief A value \c e such that analytic functions are evaluated to a tolerance of \c e.
    virtual RawFloatDP tolerance() const = 0;
    //! \brief A value \c c such that \c |a-c1| is approximately minimised.
    virtual ValueType average() const = 0;
    //! \brief A value \c c such that \c |a-c1| is approximately minimised.
    virtual NormType radius() const = 0;
    //! \brief The interval \c [c-r,c+r] where \c |a-c1|<=r.
    virtual RangeType range() const = 0;
    //! \brief An over-approximation to the norm.
    virtual NormType norm() const = 0;
};

//! \brief Interface for a unital algebra over a field with support for composition with a power series.
template<class X> class GradedAlgebraInterface
    : public virtual AlgebraInterface<X>
{
  public:
    GradedAlgebra<X> create() const;
    GradedAlgebra<X> clone() const;
  public:
    // Overrides for AlgebraInterface operations
    virtual GradedAlgebraInterface<X>* _create_copy() const = 0;
    virtual GradedAlgebraInterface<X>* _create_zero() const = 0;

    virtual Nat degree() const = 0;
    virtual const X& value() const = 0;
};

//! \brief Interface for a unital algebra over a field \a X.
template<class X> class SymbolicAlgebraInterface
    : public virtual AlgebraInterface<X>
{
  public:
    SymbolicAlgebra<X> create() const;
    SymbolicAlgebra<X> clone() const;
  public:
    // Overrides for AlgebraInterface operations
    virtual SymbolicAlgebraInterface<X>* _create_copy() const = 0;
    virtual SymbolicAlgebraInterface<X>* _create_zero() const = 0;

    using AlgebraInterface<X>::_apply;
    virtual SymbolicAlgebraInterface<X>* _apply(OperatorCode op) = 0;
};


template<class X> class VectorAlgebraInterface
    : public virtual WritableInterface
{
  public:
    virtual AlgebraInterface<X>* get_ptr(Nat i) const = 0;
    virtual Void set_ptr(Nat i, AlgebraInterface<X>*) = 0;

    Algebra<X> operator[](Nat i) const;
};


} // namespace Ariadne

#endif /* ARIADNE_ALGEBRA_INTERFACE_HPP */
