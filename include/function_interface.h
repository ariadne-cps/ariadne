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

namespace Ariadne {

typedef double Float;
class Interval;
class TaylorModel;

template<class X> class Vector;
template<class X> class Matrix;
template<class X> class Differential;

//! \brief Interface for functions whose derivatives can be computed.
class FunctionInterface {
  public:
    typedef unsigned int size_type;
    typedef unsigned short smoothness_type;

    //! \brief Virtual destructor.
    virtual ~FunctionInterface() { };
    //! \brief Create a dynamically-allocated copy.
    virtual FunctionInterface* clone() const = 0;
     
    //! \brief The smoothness of the function.
    virtual ushort smoothness() const = 0;
    //! \brief The number of arguments to the function.
    virtual uint argument_size() const = 0;
    //! \brief The number of result variables of the function.
    virtual uint result_size() const = 0;

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

    //! \brief Compute an approximation to the Jacobian derivative matrix \f$(Df)_{ij}=\partial f_i/\partial x_j\f$ of the function at the point \a x.
    virtual Matrix<Float> jacobian(const Vector<Float>& x) const = 0;
    //! \brief Compute an over-approximation to the Jacobian derivative matrix \f$(Df)_{ij}=\partial f_i/\partial x_j\f$ of the function over the domain \a x.
    virtual Matrix<Interval> jacobian(const Vector<Interval>& x) const = 0;
/*
    //! \brief Compute an approximation to all the parital derivatives \f$D^\alpha f_{i}=\partial^{|\alpha|} f_i/\partial x_\alpha\f$ of the function at the point \a x up to degree \a d.
    virtual Vector< Differential<Float> > expansion(const Vector<Float>& x, const ushort& d) const = 0;
    //! \brief Compute over-approximations to all the parital derivatives \f$D^\alpha f_{i}=\partial^{|\alpha|} f_i/\partial x_\alpha\f$ of the function over the domain \a x up to degree \a d.
    virtual Vector< Differential<Interval> > expansion(const Vector<Interval>& x, const ushort& d) const = 0;
*/
  
    //! \brief Write to an output stream.
    virtual std::ostream& write(std::ostream& os) const = 0;

    //! \brief Write to an output stream. Calls the write(std::ostream&) method to perform dynamic dispatching.
    friend std::ostream& operator<<(std::ostream& os, const FunctionInterface& f);
};

inline std::ostream& operator<<(std::ostream& os, const FunctionInterface& f) {
    return f.write(os); 
}


} // namespace Ariadne

#endif
