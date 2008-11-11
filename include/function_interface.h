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

template<class X> class Vector;
template<class X> class Matrix;
template<class X> class SparseDifferential;
template<class X> class DifferentialVector;

//! \brief Interface for functions whose derivatives can be computed.
class FunctionInterface {
  public:
    virtual ~FunctionInterface() { };
    virtual FunctionInterface* clone() const = 0;
     
    virtual ushort smoothness() const = 0;
    virtual uint argument_size() const = 0;
    virtual uint result_size() const = 0;

    virtual Vector<Float> evaluate(const Vector<Float>& x) const = 0;
    virtual Vector<Interval> evaluate(const Vector<Interval>& x) const = 0;
    virtual Matrix<Float> jacobian(const Vector<Float>& x) const = 0;
    virtual Matrix<Interval> jacobian(const Vector<Interval>& x) const = 0;
    virtual DifferentialVector< SparseDifferential<Float> > expansion(const Vector<Float>& x, const ushort& s) const = 0;
    virtual DifferentialVector< SparseDifferential<Interval> > expansion(const Vector<Interval>& x, const ushort& s) const = 0;
  
    virtual std::ostream& write(std::ostream& os) const = 0;
};

inline std::ostream& operator<<(std::ostream& os, const FunctionInterface& f) {
    return f.write(os); 
}


} // namespace Ariadne

#endif
