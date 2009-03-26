/***************************************************************************
 *            expression_interface.h
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
 
/*! \file expression_interface.h
 *  \brief Interface for real-valued expressions for which derivatives can be computed.
 */
#ifndef ARIADNE_EXPRESSION_INTERFACE_H
#define ARIADNE_EXPRESSION_INTERFACE_H

#include <iosfwd>
#include <iostream>
#include "numeric.h"

namespace Ariadne {

typedef double Float;
class Interval;
class TaylorModel;

template<class X> class Vector;
template<class X> class Matrix;
template<class X> class Differential;

//! \brief Interface for expressions whose derivatives can be computed.
class ExpressionInterface {
  public:
    //! \brief Virtual destructor.
    virtual ~ExpressionInterface() { };
    //! \brief Create a dynamically-allocated copy.
    virtual ExpressionInterface* clone() const = 0;
     
    //! \brief The smoothness of the expression.
    virtual ushort smoothness() const = 0;
    //! \brief The number of arguments to the expression.
    virtual uint argument_size() const = 0;

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

/*
    //! \brief Compute an approximation to the gradient covector \f$(Df)_{j}=\partial f/\partial x_j\f$ of the expression at the point \a x.
    virtual Vector<Float> gradient(const Vector<Float>& x) const = 0;
    //! \brief Compute an over-approximation to the Jacobian derivative matrix \f$(Df)_{ij}=\partial f/\partial x_j\f$ of the expression over the domain \a x.
    virtual Vector<Interval> gradient(const Vector<Interval>& x) const = 0;
*/
  
    //! \brief Write to an output stream.
    virtual std::ostream& write(std::ostream& os) const = 0;

    //! \brief Write to an output stream. Calls the write(std::ostream&) method to perform dynamic dispatching.
    friend std::ostream& operator<<(std::ostream& os, const ExpressionInterface& f);
};

inline std::ostream& operator<<(std::ostream& os, const ExpressionInterface& f) {
    return f.write(os); 
}


} // namespace Ariadne

#endif
