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

typedef double Float;
class Interval;
class TaylorModel;

template<class X> class Vector;
template<class X> class Matrix;
template<class X> class Differential;

class ScalarFunction;
class VectorFunction;

static const int SMOOTH=255;

//! \brief Interface for scalar functions \f$\R^n\rightarrow\R\f$.
//! See also VectorFunctionInterface.
class ScalarFunctionInterface {
  public:
    //! \brief The type used to describe the number of argument variables.
    typedef unsigned int SizeType;

    //! \brief Virtual destructor.
    virtual ~ScalarFunctionInterface() { };
    //! \brief The number of arguments to the expression.
    virtual SizeType argument_size() const = 0;

    //! \brief Compute an approximation to the value of the expression at the point \a x.
    virtual Float evaluate(const Vector<Float>& x) const = 0;
    //! \brief Compute an over-approximation to the values of the expression over the domain \a x. This method provides an <em>interval extension</em> of the expression.
    virtual Interval evaluate(const Vector<Interval>& x) const = 0;

    //! \brief Evaluate the expression over a vector of Taylor variables.
    virtual TaylorModel evaluate(const Vector<TaylorModel>& x) const = 0;
    //! \brief Evaluate the expression over a vector of differentials.
    virtual Differential<Float> evaluate(const Vector< Differential<Float> >& x) const = 0;
    //! \brief Evaluate the expression over a vector of interval differentials.
    virtual Differential<Interval> evaluate(const Vector< Differential<Interval> >& x) const = 0;

    //! \brief The derivative with respect to the \a j<sup>th</sup> coordinate.
    virtual ScalarFunction derivative(uint i) const = 0;


    //! \brief Write a brief version to an output stream.
    virtual std::ostream& repr(std::ostream& os) const = 0;
    //! \brief Write to an output stream.
    virtual std::ostream& write(std::ostream& os) const = 0;
};

//! \relates ScalarFunctionInterface
//! \brief Write to an output stream. Calls the write(std::ostream&) method to perform dynamic dispatching.
inline std::ostream& operator<<(std::ostream& os, const ScalarFunctionInterface& f) {
    return f.write(os);
}


//! \brief Interface for vector functions \f$\R^n\rightarrow\R^m\f$ whose derivatives can be computed.
//! See also ScalarFunctionInterface
class VectorFunctionInterface {
  public:
    //! \brief The type used to describe the number of argument variables.
    typedef unsigned int SizeType;

    //! \brief Virtual destructor.
    virtual ~VectorFunctionInterface() { };

    //! \brief The number of arguments to the function.
    virtual SizeType argument_size() const = 0;
    //! \brief The number of result variables of the function.
    virtual SizeType result_size() const = 0;

    //! \brief Compute an approximation to the value of the function at the point \a x.
    virtual Vector<Float> evaluate(const Vector<Float>& x) const = 0;
    //! \brief Compute an over-approximation to the values of the function over the domain \a x. This method provides an <em>interval extension</em> of the function.
    virtual Vector<Interval> evaluate(const Vector<Interval>& x) const = 0;

    //! \brief Evaluate the function over a vector of Taylor variables.
    virtual Vector<TaylorModel> evaluate(const Vector<TaylorModel>& x) const = 0;
    //! \brief Evaluate the function over a vector of differentials.
    virtual Vector< Differential<Float> > evaluate(const Vector< Differential<Float> >& x) const = 0;
    //! \brief Evaluate the function over a vector of interval differentials.
    virtual Vector< Differential<Interval> > evaluate(const Vector< Differential<Interval> >& x) const = 0;

    //! \brief Get the \a i<sup>th</sup> component function.
    virtual ScalarFunction operator[](uint i) const = 0;

    //! \brief Write to an output stream.
    virtual std::ostream& write(std::ostream& os) const = 0;
};

//! \relates VectorFunctionInterface
//! \brief Write to an output stream. Calls the write(std::ostream&) method to perform dynamic dispatching.
inline std::ostream& operator<<(std::ostream& os, const VectorFunctionInterface& f) {
    return f.write(os);
}



//! \brief An interface for scalar function models on a restricted domain.
class ScalarModelInterface {
    virtual Vector<Interval> domain() const = 0;
    virtual Interval evaluate(const Vector<Interval>&) const = 0;
};


} // namespace Ariadne

#endif
