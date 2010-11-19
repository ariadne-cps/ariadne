/***************************************************************************
 *            function_interface.h
 *
 *  Copyright 2008  Pieter Collins
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

/*! \file function_interface.h
 *  \brief Interface for functions for which derivatives can be computed.
 */
#ifndef ARIADNE_FUNCTION_INTERFACE_H
#define ARIADNE_FUNCTION_INTERFACE_H

#include <iosfwd>
#include <iostream>
#include "numeric.h"
#include "pointer.h"

namespace Ariadne {

typedef void Void;

class Float;
class Interval;
class Real;

template<class X> class Vector;
template<class X> class Matrix;
template<class X> class Differential;
template<class X> class TaylorModel;
template<class X> class Formula;

static const int SMOOTH=255;

template<class X> class ScalarFunctionInterface;
template<class X> class VectorFunctionInterface;

typedef ScalarFunctionInterface<Float> FloatScalarFunctionInterface;
typedef ScalarFunctionInterface<Interval> IntervalScalarFunctionInterface;
typedef ScalarFunctionInterface<Real> RealScalarFunctionInterface;

typedef VectorFunctionInterface<Float> FloatVectorFunctionInterface;
typedef VectorFunctionInterface<Interval> IntervalVectorFunctionInterface;
typedef VectorFunctionInterface<Real> RealVectorFunctionInterface;

template<class X> class ScalarFunction;
typedef ScalarFunction<Float> FloatScalarFunction;
typedef ScalarFunction<Interval> IntervalScalarFunction;
typedef ScalarFunction<Real> RealScalarFunction;

template<class X> class VectorFunction;
typedef VectorFunction<Float> FloatVectorFunction;
typedef VectorFunction<Interval> IntervalVectorFunction;
typedef VectorFunction<Real> RealVectorFunction;

template<>
class ScalarFunctionInterface<Void>
{
  public:
    //! \brief The type used to describe the number of argument variables.
    typedef unsigned int SizeType;

    //! \brief Virtual destructor.
    virtual ~ScalarFunctionInterface() { };
    //! \brief The number of arguments to the expression.
    virtual SizeType argument_size() const = 0;

    //! \brief Write a brief version to an output stream.
    virtual std::ostream& repr(std::ostream& os) const = 0;
    //! \brief Write to an output stream.
    virtual std::ostream& write(std::ostream& os) const = 0;
};

//! \ingroup FunctionModule
//! \brief Interface for scalar functions \f$\mathbb{F}^n\rightarrow\mathbb{F}\f$ which can only be evaluated approximately.
//! \sa \ref VectorFunctionInterface.
template<>
class ScalarFunctionInterface<Float>
    : public ScalarFunctionInterface<Void>
{
  public:
    //! \brief Compute an approximation to the value of the function at the point \a x.
    virtual Float evaluate(const Vector<Float>& x) const = 0;
    //! \brief Evaluate the function over a vector of differentials.
    virtual Differential<Float> evaluate(const Vector< Differential<Float> >& x) const = 0;
    //! \brief Evaluate the function over a vector of approximate Taylor models.
    virtual TaylorModel<Float> evaluate(const Vector< TaylorModel<Float> >& x) const = 0;
};

//! \ingroup FunctionModule
//! \brief Interface for scalar functions \f$\mathbb{I}^n\rightarrow\mathbb{I}\f$ which can be evaluated over intervals.
//! \sa \ref VectorFunctionInterface.
template<>
class ScalarFunctionInterface<Interval>
    : public ScalarFunctionInterface<Float>
{
  public:
    //! \brief The type used to describe the number of argument variables.
    typedef unsigned int SizeType;

    using ScalarFunctionInterface<Float>::evaluate;

    //! \brief Compute an over-approximation to the values of the function over the domain \a x. This method provides an <em>interval extension</em> of the function.
    virtual Interval evaluate(const Vector<Interval>& x) const = 0;
    //! \brief Evaluate the function over a vector of interval differentials.
    virtual Differential<Interval> evaluate(const Vector< Differential<Interval> >& x) const = 0;
    //! \brief Evaluate the function over a vector of Taylor models with interval error.
    virtual TaylorModel<Interval> evaluate(const Vector< TaylorModel<Interval> >& x) const = 0;

    //! \brief Apply the function to a formula. Can be used to obtain a tree structure from the function.
    virtual Formula<Interval> evaluate(const Vector< Formula<Interval> >& x) const = 0;

};

//! \ingroup FunctionModule
//! \brief Interface for scalar functions \f$\R^n\rightarrow\R\f$ which can be evaluated exactly.
//! \sa \ref VectorFunctionInterface.
template<>
class ScalarFunctionInterface<Real>
    : public ScalarFunctionInterface<Interval>
{
  public:
    using ScalarFunctionInterface<Interval>::evaluate;

    //! \brief Evaluate over computable reals.
    virtual Real evaluate(const Vector<Real>& x) const = 0;

    //! \brief The derivative with respect to the \a j<sup>th</sup> coordinate.
    virtual ScalarFunction<Real> derivative(uint i) const = 0;
};

//! \relates ScalarFunctionInterface
//! \brief Write to an output stream. Calls the write(std::ostream&) method to perform dynamic dispatching.
inline std::ostream& operator<<(std::ostream& os, const ScalarFunctionInterface<Void>& f) {
    return f.write(os);
}


//! \ingroup FunctionModule
//! \brief Interface for vector functions \f$\F^n\rightarrow\F^m\f$ whose derivatives can be computed.
//! \sa \ref ScalarFunctionInterface
template<>
class VectorFunctionInterface<Void>
{
  public:
    //! \brief The type used to describe the number of argument variables.
    typedef unsigned int SizeType;

    //! \brief Virtual destructor.
    virtual ~VectorFunctionInterface() { };

    //! \brief The number of arguments to the function.
    virtual SizeType argument_size() const = 0;
    //! \brief The number of result variables of the function.
    virtual SizeType result_size() const = 0;

    //! \brief Write to an output stream.
    virtual std::ostream& write(std::ostream& os) const = 0;
};
//! \ingroup FunctionModule
//! \brief Interface for vector functions \f$\F^n\rightarrow\F^m\f$ whose derivatives can be computed.
//! \sa \ref ScalarFunctionInterface
template<>
class VectorFunctionInterface<Float>
    : public VectorFunctionInterface<Void>
{
  public:
    //! \brief Compute an approximation to the value of the function at the point \a x.
    virtual Vector<Float> evaluate(const Vector<Float>& x) const = 0;
    //! \brief Evaluate the function over a vector of differentials.
    virtual Vector< Differential<Float> > evaluate(const Vector< Differential<Float> >& x) const = 0;
    //! \brief Evaluate the function over a vector of approximate Taylor models.
    virtual Vector< TaylorModel<Float> > evaluate(const Vector< TaylorModel<Float> >& x) const = 0;
};

//! \ingroup FunctionModule
//! \brief Interface for vector functions \f$\I^n\rightarrow\I^m\f$ whose derivatives can be computed.
//! \sa \ref ScalarFunctionInterface
template<>
class VectorFunctionInterface<Interval>
    : public VectorFunctionInterface<Float>
{
  public:
    using VectorFunctionInterface<Float>::evaluate;

    //! \brief Compute an over-approximation to the values of the function over the domain \a x. This method provides an <em>interval extension</em> of the function.
    virtual Vector<Interval> evaluate(const Vector<Interval>& x) const = 0;
    //! \brief Evaluate the function over a vector of interval differentials.
    virtual Vector< Differential<Interval> > evaluate(const Vector< Differential<Interval> >& x) const = 0;
    //! \brief Evaluate the function over a vector of Taylor models with interval error.
    virtual Vector< TaylorModel<Interval> > evaluate(const Vector< TaylorModel<Interval> >& x) const = 0;

    //! \brief Evaluate the function over a vector of formulae.
    virtual Vector< Formula<Interval> > evaluate(const Vector< Formula<Interval> >& x) const = 0;
};

//! \ingroup FunctionModule
//! \brief Interface for vector functions \f$\R^n\rightarrow\R^m\f$ whose derivatives can be computed.
//! \sa \ref ScalarFunctionInterface
template<>
class VectorFunctionInterface<Real>
    : public VectorFunctionInterface<Interval>
{
  public:
    using VectorFunctionInterface<Interval>::evaluate;

    //! \brief Evaluate over computable reals.
    virtual Vector<Real> evaluate(const Vector<Real>& x) const = 0;
    //! \brief Get the \a i<sup>th</sup> component function.
    virtual ScalarFunction<Real> operator[](uint i) const = 0;
};

//! \relates VectorFunctionInterface
//! \brief Write to an output stream. Calls the write(std::ostream&) method to perform dynamic dispatching.
inline std::ostream& operator<<(std::ostream& os, const VectorFunctionInterface<Void>& f) {
    return f.write(os);
}



//! \brief An interface for scalar function models on a restricted domain.
class ScalarModelInterface {
    virtual Vector<Interval> domain() const = 0;
    virtual Interval evaluate(const Vector<Interval>&) const = 0;
};


} // namespace Ariadne

#endif
