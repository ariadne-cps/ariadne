/***************************************************************************
 *            constraint.h
 *
 *  Copyright 2009  Pieter Collins
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

/*! \file constraint.h
 *  \brief Linear and nonlinear constraints.
 */

#ifndef ARIADNE_CONSTRAINT_H
#define ARIADNE_CONSTRAINT_H

#include "numeric.h"
#include "function.h"
#include "taylor_function.h"

namespace Ariadne {

class NonlinearConstraint {
  public:
    NonlinearConstraint(ScalarFunction const& f, Interval const& b) : _function(f), _bounds(b) { }
    ScalarFunction const& function() const { return this->_function; }
    Interval const& bounds() const { return this->_bounds; }
  private:
    ScalarFunction _function;
    Interval _bounds;
};

inline NonlinearConstraint operator<=(const Float& c, const ScalarFunction& f) {
    return NonlinearConstraint(f,Interval(c,inf<Float>()));
}

inline NonlinearConstraint operator<=(const ScalarFunction& f, const Float& c) {
    return NonlinearConstraint(f,Interval(-inf<Float>(),c));
}

inline NonlinearConstraint operator>=(const ScalarFunction& f, const Float& c) {
    return NonlinearConstraint(f,Interval(c,+inf<Float>()));
}

inline NonlinearConstraint operator==(const ScalarFunction& f, double c) {
    return NonlinearConstraint(f,Interval(c));
}

inline NonlinearConstraint operator==(const ScalarFunction& f, const Real& c) {
    return NonlinearConstraint(f,Interval(c));
}

inline NonlinearConstraint operator==(const ScalarFunction& f, const Float& c) {
    return NonlinearConstraint(f,Interval(c));
}

inline NonlinearConstraint operator==(const ScalarFunction& f, const Interval& c) {
    return NonlinearConstraint(f,Interval(c));
}

inline NonlinearConstraint operator<=(const ScalarTaylorFunction& f, const Float& c) {
    return NonlinearConstraint(f.function(),Interval(-inf<Float>(),c));
}

inline NonlinearConstraint operator>=(const ScalarTaylorFunction& f, const Float& c) {
    return NonlinearConstraint(f.function(),Interval(c,+inf<Float>()));
}

inline NonlinearConstraint operator<=(const ScalarTaylorFunction& f, const ScalarTaylorFunction& g) {
    return NonlinearConstraint((f-g).function(),Interval(-inf<Float>(),0.0));
}

inline NonlinearConstraint operator>=(const ScalarTaylorFunction& f, const ScalarTaylorFunction& g) {
    return NonlinearConstraint((f-g).function(),Interval(0.0,inf<Float>()));
}

inline NonlinearConstraint operator<=(const NonlinearConstraint& nc, const Float& c) {
    return NonlinearConstraint(nc.function(),intersection(nc.bounds(),Interval(-inf<Float>(),c)));
}

std::ostream& operator<<(std::ostream& os, const NonlinearConstraint& c);


} //namespace Ariadne

#endif
