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

namespace Ariadne {

template<class X, class R>
class NonlinearConstraint {
    typedef ScalarFunction<X> FunctionType;
    typedef R BoundType;
  public:
    NonlinearConstraint(R const& l, ScalarFunction<X> const& f, R const& u) : _function(f), _lower_bound(l), _upper_bound(u) { ARIADNE_ASSERT(l<=u); }
    NonlinearConstraint(ScalarFunction<X> const& f, R const& x) : _function(f), _lower_bound(x), _upper_bound(x) { }
    template<class XX,class RR> NonlinearConstraint(const NonlinearConstraint<XX,RR>& c)
        : _function(static_cast<ScalarFunction<X> >(c.function())), _lower_bound(c.lower_bound()), _upper_bound(c.upper_bound()) { }
    ScalarFunction<X> const& function() const { return this->_function; }
    R const& lower_bound() const { return this->_lower_bound; }
    R const& upper_bound() const { return this->_upper_bound; }
    const Interval bounds() const { return Interval(this->_lower_bound,this->_upper_bound); }
  private:
    ScalarFunction<X> _function;
    R _lower_bound;
    R _upper_bound;
};

typedef NonlinearConstraint<Real,Real> RealNonlinearConstraint;
typedef NonlinearConstraint<Interval,Float> IntervalNonlinearConstraint;

inline RealNonlinearConstraint operator<=(const Real& c, const RealScalarFunction& f) {
    return RealNonlinearConstraint(c,f,infinity);
}

inline RealNonlinearConstraint operator>=(const Real& c, const RealScalarFunction& f) {
    return RealNonlinearConstraint(-infinity,f,c);
}

inline RealNonlinearConstraint operator<=(const RealScalarFunction& f, const Real& c) {
    return RealNonlinearConstraint(-infinity,f,c);
}

inline RealNonlinearConstraint operator>=(const RealScalarFunction& f, const Real& c) {
    return RealNonlinearConstraint(c,f,infinity);
}

inline RealNonlinearConstraint operator==(const RealScalarFunction& f, const Real& c) {
    return RealNonlinearConstraint(f,c);
}

inline RealNonlinearConstraint operator<=(const RealScalarFunction& f, double c) {
    return RealNonlinearConstraint(-infinity,f,Real(c));
}

inline RealNonlinearConstraint operator>=(const RealScalarFunction& f, double c) {
    return RealNonlinearConstraint(Real(c),f,infinity);
}

inline RealNonlinearConstraint operator==(const RealScalarFunction& f, double c) {
    return RealNonlinearConstraint(f,Real(c));
}


inline RealNonlinearConstraint operator<=(const RealNonlinearConstraint& nc, const Real& c) {
    ARIADNE_ASSERT(Float(nc.upper_bound())==inf<Float>());
    return RealNonlinearConstraint(nc.lower_bound(),nc.function(),c);
}


inline IntervalNonlinearConstraint operator<=(const IntervalScalarFunction& f, const Float& c) {
    return IntervalNonlinearConstraint(-inf<Float>(),f,c);
}

inline IntervalNonlinearConstraint operator>=(const IntervalScalarFunction& f, const Float& c) {
    return IntervalNonlinearConstraint(c,f,+inf<Float>());
}

inline IntervalNonlinearConstraint operator==(const IntervalScalarFunction& f, const Float& c) {
    return IntervalNonlinearConstraint(c,f,c);
}

inline IntervalNonlinearConstraint operator<=(const IntervalScalarFunction& f, double c) {
    return IntervalNonlinearConstraint(-inf<Float>(),f,Float(c));
}

inline IntervalNonlinearConstraint operator>=(const IntervalScalarFunction& f, double c) {
    return IntervalNonlinearConstraint(-inf<Float>(),f,c);
}

inline IntervalNonlinearConstraint operator==(const IntervalScalarFunction& f, double c) {
    return IntervalNonlinearConstraint(Float(c),f,Float(c));
}

inline IntervalNonlinearConstraint operator<=(const IntervalScalarFunction& f1, const IntervalScalarFunction& f2) {
    return (f1-f2) <= 0.0;
}

inline IntervalNonlinearConstraint operator>=(const IntervalScalarFunction& f1, const IntervalScalarFunction& f2) {
    return (f1-f2) >= 0.0;
}


template<class X, class R> std::ostream& operator<<(std::ostream& os, const NonlinearConstraint<X,R>& c) {
    return os << c.lower_bound() << "<=" << c.function() << "<=" << c.upper_bound();
}


} //namespace Ariadne

#endif
