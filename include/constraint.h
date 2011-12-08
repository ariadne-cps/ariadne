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


template<class F, class R>
class Constraint {
    typedef F FunctionType;
    typedef R BoundType;
//    typedef typename IntervalOfType<Real>::Type IntervalBoundsType;
  public:
    Constraint(BoundType const& l, FunctionType const& f, BoundType const& u) : _function(f), _lower_bound(l), _upper_bound(u) { ARIADNE_ASSERT(l<=u); }
    Constraint(FunctionType const& f, BoundType const& x) : _function(f), _lower_bound(x), _upper_bound(x) { }
    template<class FF,class RR> Constraint(const Constraint<FF,RR>& c)
        : _function(static_cast<F>(c.function())), _lower_bound(c.lower_bound()), _upper_bound(c.upper_bound()) { }
    void set_function(const FunctionType& f) { this->_function = f; }
    FunctionType& function() { return this->_function; }
    FunctionType const& function() const { return this->_function; }
    Nat argument_size() const { return this->_function.argument_size(); }
    BoundType const& lower_bound() const { return this->_lower_bound; }
    BoundType const& upper_bound() const { return this->_upper_bound; }
    const Interval bounds() const { return Interval(this->_lower_bound,this->_upper_bound); }
  private:
    F _function;
    R _lower_bound;
    R _upper_bound;
};

typedef Constraint<RealScalarFunction,Real> RealConstraint;
typedef Constraint<IntervalScalarFunction,Float> IntervalConstraint;

inline RealConstraint operator<=(const Real& c, const RealScalarFunction& f) {
    return RealConstraint(c,f,infinity);
}

inline RealConstraint operator>=(const Real& c, const RealScalarFunction& f) {
    return RealConstraint(-infinity,f,c);
}

inline RealConstraint operator<=(const RealScalarFunction& f, const Real& c) {
    return RealConstraint(-infinity,f,c);
}

inline RealConstraint operator>=(const RealScalarFunction& f, const Real& c) {
    return RealConstraint(c,f,infinity);
}

inline RealConstraint operator==(const RealScalarFunction& f, const Real& c) {
    return RealConstraint(f,c);
}

inline RealConstraint operator<=(const RealScalarFunction& f, double c) {
    return RealConstraint(-infinity,f,Real(c));
}

inline RealConstraint operator>=(const RealScalarFunction& f, double c) {
    return RealConstraint(Real(c),f,infinity);
}

inline RealConstraint operator==(const RealScalarFunction& f, double c) {
    return RealConstraint(f,Real(c));
}


inline RealConstraint operator<=(const RealConstraint& nc, const Real& c) {
    ARIADNE_ASSERT(Float(nc.upper_bound())==inf);
    return RealConstraint(nc.lower_bound(),nc.function(),c);
}


inline IntervalConstraint operator<=(const Float& c, const IntervalScalarFunction& f) {
    return IntervalConstraint(c,f,+inf);
}

inline IntervalConstraint operator<=(const IntervalScalarFunction& f, const Float& c) {
    return IntervalConstraint(-inf,f,c);
}

inline IntervalConstraint operator>=(const IntervalScalarFunction& f, const Float& c) {
    return IntervalConstraint(c,f,+inf);
}

inline IntervalConstraint operator==(const IntervalScalarFunction& f, const Float& c) {
    return IntervalConstraint(c,f,c);
}

inline IntervalConstraint operator<=(const IntervalScalarFunction& f, double c) {
    return IntervalConstraint(-inf,f,Float(c));
}

inline IntervalConstraint operator>=(const IntervalScalarFunction& f, double c) {
    return IntervalConstraint(-inf,f,c);
}

inline IntervalConstraint operator==(const IntervalScalarFunction& f, double c) {
    return IntervalConstraint(Float(c),f,Float(c));
}

inline IntervalConstraint operator<=(const IntervalScalarFunction& f1, const IntervalScalarFunction& f2) {
    return (f1-f2) <= 0.0;
}

inline IntervalConstraint operator>=(const IntervalScalarFunction& f1, const IntervalScalarFunction& f2) {
    return (f1-f2) >= 0.0;
}

inline IntervalConstraint operator<=(const IntervalConstraint& nc, const Float& c) {
    ARIADNE_ASSERT(nc.upper_bound()==inf);
    return IntervalConstraint(nc.lower_bound(),nc.function(),c);
}


template<class X, class R> std::ostream& operator<<(std::ostream& os, const Constraint<X,R>& c) {
    return os << c.lower_bound() << "<=" << c.function() << "<=" << c.upper_bound();
}


} //namespace Ariadne

#endif
