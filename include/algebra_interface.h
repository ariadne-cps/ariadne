/***************************************************************************
 *            algebra_interface.h
 *
 *  Copyright 2010  Pieter Collins
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
#include "numeric.h"
#include "pointer.h"
#include "operators.h"

namespace Ariadne {

class Float;
class Interval;
class Real;

template<class X> class Vector;
template<class X> class Matrix;
template<class X> class Differential;

template<class X> class Series;

template<class X> class AlgebraInterface;
template<class X> class Algebra;
template<class X> class NormedAlgebra;
template<class X> class GradedAlgebra;
template<class X> class SymbolicAlgebra;

class WritableInterface {
    friend std::ostream& operator<<(std::ostream& os, const WritableInterface& w);
  public:
    virtual std::ostream& write(std::ostream&) const = 0;
};
inline std::ostream& operator<<(std::ostream& os, const WritableInterface& w) {
    w.write(os); return os; }

struct Ball {
    explicit Ball(Float r) : _radius(r) { }
    Float _radius;
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
    //! \brief Create a dynamically-allocated copy.
    virtual AlgebraInterface<X>* _clone() const = 0;
    //! \brief Create the zero element in the same algebra as the current object.
    virtual AlgebraInterface<X>* _create() const = 0;

    //! \brief Add a constant numerical scalar \c r+=c .
    virtual void _iadd(const X& c) = 0;
    //! \brief Multiply by a numerical scalar \c r*=c .
    virtual void _imul(const X& c) = 0;
    //! \brief Scalar multiply and add \c r+=c*x .
    virtual void _isma(const X& c, const AlgebraInterface<X>& x) = 0;
    //! \brief Fused multiply and add \c r+=x1*x2 .
    virtual void _ifma(const AlgebraInterface<X>& x1, const AlgebraInterface<X>& x2) = 0;
};

//! \brief Interface for a normed unital algebra over a field \a X.
template<class X> class NormedAlgebraInterface
    : public virtual AlgebraInterface<X>
{
  public:
    NormedAlgebra<X> create() const;
    NormedAlgebra<X> clone() const;
  public:
    // Overrides for AlgebraInterface operations
    virtual NormedAlgebraInterface<X>* _clone() const = 0;
    virtual NormedAlgebraInterface<X>* _create() const = 0;
    virtual NormedAlgebraInterface<X>* _create_ball(Float r) const = 0;

    //! \brief A value \c e such that analytic functions are evaluated to a tolerance of \c e.
    virtual Float tolerance() const = 0;
    //! \brief A value \c c such that \c |a-c1| is approximately minimised.
    virtual Float average() const = 0;
    //! \brief A value \c r such that \c |a-c1|<=r.
    virtual Float radius() const = 0;
    //! \brief An over-approximation to the norm.
    virtual Float norm() const = 0;
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

    virtual uint degree() const = 0;
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

    virtual SymbolicAlgebraInterface<X>* _apply(Operator op) = 0;
};


template<class X> class VectorAlgebraInterface
    : public virtual WritableInterface
{
  public:
    virtual AlgebraInterface<X>* get_ptr(uint i) const = 0;
    virtual void set_ptr(uint i, AlgebraInterface<X>*) = 0;

    Algebra<X> operator[](uint i) const;
};


} // namespace Ariadne

#endif /* ARIADNE_ALGEBRA_INTERFACE_H */
