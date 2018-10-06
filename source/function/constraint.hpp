/***************************************************************************
 *            constraint.hpp
 *
 *  Copyright 2009-17  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file constraint.hpp
 *  \brief Linear and nonlinear constraints.
 */

#ifndef ARIADNE_CONSTRAINT_HPP
#define ARIADNE_CONSTRAINT_HPP

#include "../numeric/numeric.hpp"
#include "../function/function.hpp"
#include "../utility/metaprogramming.hpp"

namespace Ariadne {

template<class F, class R>
class Constraint {
  public:
    typedef F FunctionType;
    typedef R BoundType;
    typedef R UpperBoundType;
    typedef decltype(-std::declval<R>()) LowerBoundType;
//    typedef typename IntervalOfType<Real>::Type IntervalBoundsType;
  public:
    Constraint(LowerBoundType const& l, FunctionType const& f, UpperBoundType const& u)
        : _function(f), _lower_bound(l), _upper_bound(u) { ARIADNE_ASSERT_MSG(decide(l<=u),"f="<<f<<"\nl="<<l<<", u="<<u); }

    Constraint(FunctionType const& f, BoundType const& x)
        : _function(f), _lower_bound(x), _upper_bound(x) { }

    // FIXME: Should require convertibility of RR to R
    template< class FF, class RR, EnableIf<IsConvertible<FF,F>> =dummy, EnableIf<IsConstructible<R,RR>> =dummy>
    Constraint(const Constraint<FF,RR>& c)
        : _function(static_cast<F>(c.function())), _lower_bound(c.lower_bound()), _upper_bound(c.upper_bound()) { }

    template< class FF, class RR, EnableIf<IsConvertible<FF,F>> =dummy, EnableIf<IsConstructible<R,RR,DoublePrecision>> =dummy>
    Constraint(const RR& l, const FF& f, const RR& u)
        : _function(static_cast<F>(f)), _lower_bound(l,dp), _upper_bound(u,dp) { }

    template< class FF, class RR, EnableIf<IsConvertible<FF,F>> =dummy, EnableIf<IsConstructible<R,RR,DoublePrecision>> =dummy>
    Constraint(const Constraint<FF,RR>& c)
        : _function(static_cast<F>(c.function())), _lower_bound(c.lower_bound(),dp), _upper_bound(c.upper_bound(),dp) { }

    Void set_function(const FunctionType& f) { this->_function = f; }
    FunctionType& function() { return this->_function; }
    FunctionType const& function() const { return this->_function; }
    Nat argument_size() const { return this->_function.argument_size(); }
    LowerBoundType const& lower_bound() const { return this->_lower_bound; }
    UpperBoundType const& upper_bound() const { return this->_upper_bound; }

    // FIXME: This function should not be used as it breaks type safety
    const Interval<FloatDPValue> bounds() const { DoublePrecision pr; return cast_exact(Interval<FloatDPUpperBound>({_lower_bound,pr},{_upper_bound,pr})); }
  private:
    F _function;
    R _lower_bound;
    R _upper_bound;
};

typedef Constraint<RealScalarFunction,Real> RealConstraint;
typedef Constraint<EffectiveScalarFunction,EffectiveNumber> EffectiveConstraint;
typedef Constraint<ValidatedScalarFunction,ValidatedNumber> ValidatedConstraint;
typedef Constraint<ValidatedScalarFunction,ExactNumber> ValidatedExactConstraint;

template<class X, class R> OutputStream& operator<<(OutputStream& os, const Constraint<X,R>& c) {
    return os << c.lower_bound() << "<=" << c.function() << "<=" << c.upper_bound();
}

inline EffectiveConstraint operator<=(const EffectiveNumber& c, const EffectiveScalarFunction& f) {
    return EffectiveConstraint(c,f,infinity);
}

inline EffectiveConstraint operator>=(const EffectiveNumber& c, const EffectiveScalarFunction& f) {
    return EffectiveConstraint(-infinity,f,c);
}

inline EffectiveConstraint operator<=(const EffectiveScalarFunction& f, const EffectiveNumber& c) {
    return EffectiveConstraint(-infinity,f,c);
}

inline EffectiveConstraint operator>=(const EffectiveScalarFunction& f, const EffectiveNumber& c) {
    return EffectiveConstraint(c,f,infinity);
}

inline EffectiveConstraint operator==(const EffectiveScalarFunction& f, const EffectiveNumber& c) {
    return EffectiveConstraint(f,c);
}

inline EffectiveConstraint operator<=(const EffectiveScalarFunction& f, double c) {
    return f <= Dyadic(c);
}

inline EffectiveConstraint operator>=(const EffectiveScalarFunction& f, double c) {
    return f >= Dyadic(c);
}

inline EffectiveConstraint operator==(const EffectiveScalarFunction& f, double c) {
    return f == Dyadic(c);
}

inline EffectiveConstraint operator<=(const EffectiveConstraint& nc, const EffectiveNumber& c) {
    ARIADNE_ASSERT(decide(nc.upper_bound()==infty));
    return EffectiveConstraint(nc.lower_bound(),nc.function(),c);
}


inline ValidatedExactConstraint operator<=(const ExactNumber& c, const ValidatedScalarFunction& f) {
    return ValidatedExactConstraint(c,f,ExactNumber(+infty));
}

inline ValidatedExactConstraint operator<=(const ValidatedScalarFunction& f, const ExactNumber& c) {
    return ValidatedExactConstraint(-infty,f,c);
}

inline ValidatedExactConstraint operator>=(const ValidatedScalarFunction& f, const ExactNumber& c) {
    return ValidatedExactConstraint(c,f,+infty);
}

inline ValidatedExactConstraint operator==(const ValidatedScalarFunction& f, const ExactNumber& c) {
    return ValidatedExactConstraint(c,f,c);
}

inline ValidatedExactConstraint operator<=(const ValidatedExactConstraint& nc, const ExactNumber& c) {
    ARIADNE_ASSERT(nc.upper_bound()==infty);
    return ValidatedExactConstraint(nc.lower_bound(),nc.function(),c);
}


inline ValidatedExactConstraint operator<=(const ValidatedScalarFunction& f1, const ValidatedScalarFunction& f2) {
    return (f1-f2) <= ExactNumber(0);
}

inline ValidatedExactConstraint operator>=(const ValidatedScalarFunction& f1, const ValidatedScalarFunction& f2) {
    return (f1-f2) >= ExactNumber(0);
}


inline ValidatedConstraint operator<=(const ValidatedNumber& c, const ValidatedScalarFunction& f) {
    return ValidatedConstraint(c,f,ExactNumber(+infty));
}

inline ValidatedConstraint operator<=(const ValidatedConstraint& nc, const ValidatedNumber& c) {
    ARIADNE_ASSERT(decide(nc.upper_bound()==infty));
   return ValidatedConstraint(nc.lower_bound(),nc.function(),c);
}



} //namespace Ariadne

#endif
