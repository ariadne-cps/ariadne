/***************************************************************************
 *            algebra_interface.h
 *
 *  Copyright 2010-15  Pieter Collins
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

/*! \file algebra_interface.h
 *  \brief Interface for function algebras.
 */

#ifndef ARIADNE_ALGEBRA_INTERFACE_H
#define ARIADNE_ALGEBRA_INTERFACE_H

#include <iosfwd>
#include <iostream>

#include "utility/writable.h"
#include "numeric/numeric.h"
#include "utility/pointer.h"
#include "expression/operators.h"

namespace Ariadne {

template<class X> class AlgebraInterface;
template<class X> class Algebra;
template<class X> class NormedAlgebra;
template<class X> class GradedAlgebra;
template<class X> class SymbolicAlgebra;

template<class X> struct AlgebraTraits;

template<> struct AlgebraTraits<ApproximateFloat> {
    typedef ApproximateFloat ValueType;
    typedef ApproximateInterval RangeType;
    typedef ApproximateFloat NormType;
    typedef ApproximateNumber NumericType;
};

template<> struct AlgebraTraits<ValidatedFloat> {
    typedef ExactFloat ValueType;
    typedef UpperInterval RangeType;
    typedef ErrorFloat NormType;
    typedef ValidatedNumber NumericType;
};

template<> struct AlgebraTraits<Real> {
    typedef ExactFloat ValueType;
    typedef UpperInterval RangeType;
    typedef ErrorFloat NormType;
    typedef Real NumericType;
};


struct Ball {
    explicit Ball(ErrorType r) : _radius(r) { }
    ErrorType _radius;
};

//! \brief Interface for a unital algebra over a field \a X.
template<class X> class AlgebraInterface
    : public virtual WritableInterface
{
    friend class Algebra<X>;
  public:
    typedef X NumericType;
  public:
    Algebra<X> create() const;
    Algebra<X> clone() const;
  public:
    //! \brief Virtual destructor.
    virtual ~AlgebraInterface<X>() { }
    //! \brief Create a dynamically-allocated copy.
    virtual AlgebraInterface<X>* _clone() const = 0;
    //! \brief Create the zero element in the same algebra as the current object.
    virtual AlgebraInterface<X>* _create() const = 0;

    //! \brief Add a constant numerical scalar \c r+=c .
    virtual Void _iadd(const X& c) = 0;
    //! \brief Multiply by a numerical scalar \c r*=c .
    virtual Void _imul(const X& c) = 0;
    //! \brief Scalar multiply and add \c r+=c*x .
    virtual Void _isma(const X& c, const AlgebraInterface<X>& x) = 0;
    //! \brief Fused multiply and add \c r+=x1*x2 .
    virtual Void _ifma(const AlgebraInterface<X>& x1, const AlgebraInterface<X>& x2) = 0;
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
    virtual NormedAlgebraInterface<X>* _clone() const = 0;
    virtual NormedAlgebraInterface<X>* _create() const = 0;
    virtual NormedAlgebraInterface<X>* _create_constant(X c) const = 0;
    virtual NormedAlgebraInterface<X>* _create_ball(ErrorType r) const = 0;

    //! \brief A value \c e such that analytic functions are evaluated to a tolerance of \c e.
    virtual RawFloat tolerance() const = 0;
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
    virtual GradedAlgebraInterface<X>* _clone() const = 0;
    virtual GradedAlgebraInterface<X>* _create() const = 0;

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
    virtual SymbolicAlgebraInterface<X>* _clone() const = 0;
    virtual SymbolicAlgebraInterface<X>* _create() const = 0;

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

#endif /* ARIADNE_ALGEBRA_INTERFACE_H */
