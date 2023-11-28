/***************************************************************************
 *            function/constraint.hpp
 *
 *  Copyright  2009-20  Pieter Collins
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

/*! \file function/constraint.hpp
 *  \brief Linear and nonlinear constraints.
 */

#ifndef ARIADNE_CONSTRAINT_HPP
#define ARIADNE_CONSTRAINT_HPP

#include "numeric/numeric.hpp"
#include "function/function.hpp"
#include "utility/metaprogramming.hpp"

namespace Ariadne {

template<class F, class R> class Constraint;
template<class F, class R> class LowerConstraint;
template<class F, class R> class UpperConstraint;

//! \ingroup GeometryModule
//! \brief A constraint defined by requiring values of a function \f$f\f$ to lie in a range \f$[l:u]\f$.
//! i.e. restricts point \f$x\f$ to satsify \f$f(x)\in [l:u]\f$.
//!    \param F The type of the function
//!    \param U The type of the upper bound of the range of values.
template<class F, class U>
class Constraint {
  public:
    typedef F FunctionType; //!< <p/>
    typedef U BoundType; //!< <p/>
    typedef U UpperBoundType; //!< <p/>
    typedef NegationType<U> LowerBoundType; //!< <p/>
//    typedef typename IntervalOfType<Real>::Type IntervalBoundsType;
  public:
    //! <p/>
    Constraint(LowerBoundType const& l, FunctionType const& f, UpperBoundType const& u)
        : _function(f), _lower_bound(l), _upper_bound(u) { ARIADNE_ASSERT_MSG(decide(l<=u),"f="<<f<<"\nl="<<l<<", u="<<u); }

    //! <p/>
    Constraint(FunctionType const& f, BoundType const& x)
        : _function(f), _lower_bound(x), _upper_bound(x) { }

    template<ConvertibleTo<F> FF, ConvertibleTo<U> UU> Constraint(const Constraint<FF,UU>& c)
        : _function(static_cast<F>(c.function())), _lower_bound(c.lower_bound()), _upper_bound(c.upper_bound()) { }

    Void set_function(const FunctionType& f) { this->_function = f; }
    FunctionType& function() { return this->_function; }
    //! <p/>
    FunctionType const& function() const { return this->_function; }
    //! <p/>
    SizeType argument_size() const { return this->_function.argument_size(); }
    //! <p/>
    LowerBoundType const& lower_bound() const { return this->_lower_bound; }
    //! <p/>
    UpperBoundType const& upper_bound() const { return this->_upper_bound; }

    //! <p/>
    const Interval<U> bounds() const {
        return Interval<U>(this->_lower_bound,this->_upper_bound); }

    // FIXME: This function should not be used as it breaks type safety
    const Interval<FloatDP> cast_exact_bounds() const {
        if constexpr (Constructible<FloatDPUpperBound,U,DoublePrecision>) {
            return cast_exact(Interval<FloatDPUpperBound>({_lower_bound,dp},{_upper_bound,dp}));
        } else {
            return cast_exact(Interval<FloatDPApproximation>({_lower_bound,dp},{_upper_bound,dp}));
        }
    }

  private:
    F _function;
    U _lower_bound;
    U _upper_bound;
};

template<class F, class U> OutputStream& operator<<(OutputStream& os, const Constraint<F,U>& c) {
    return os << c.lower_bound() << "<=" << c.function() << "<=" << c.upper_bound();
}


template<class F, class U> class LowerConstraint : public Constraint<F,U> {
    F _function; U _lower_bound;
  public:
    LowerConstraint(U const& l, F const& f) : Constraint<F,U>(l,f,+infty) { }
    friend OutputStream& operator<<(OutputStream& os, const LowerConstraint<F,U>& c) {
        return os << c.function() << ">=" << c.lower_bound(); }
};

template<class F, class U> class UpperConstraint : public Constraint<F,U> {
    F _function; U _upper_bound;
  public:
    UpperConstraint(F const& f, U const& u) : Constraint<F,U>(-infty,f,u) { }
    friend OutputStream& operator<<(OutputStream& os, const UpperConstraint<F,U>& c) {
        return os << c.function() << "<=" << c.upper_bound(); }
};

template<class F, class U> class EqualityConstraint : public Constraint<F,U> {
    F _function; U _value;
  public:
    EqualityConstraint(F const& f, U const& v) : Constraint<F,U>(v,f,v) { }
    friend OutputStream& operator<<(OutputStream& os, const EqualityConstraint<F,U>& c) {
        return os << c.function() << "==" << c._value; }
};


//!@{
//! \relates Constraint
//! \name Type synonyms
template<class P> using ConstraintType = Constraint<ScalarMultivariateFunction<P>,Number<P>>; //!< <p/>
template<class P> using LowerConstraintType = LowerConstraint<ScalarMultivariateFunction<P>,Number<P>>; //!< <p/>
template<class P> using UpperConstraintType = UpperConstraint<ScalarMultivariateFunction<P>,Number<P>>; //!< <p/>
template<class P> using EqualityConstraintType = EqualityConstraint<ScalarMultivariateFunction<P>,Number<P>>; //!< <p/>

using RealConstraint = Constraint<RealScalarMultivariateFunction,Real>; //!< <p/>
using EffectiveConstraint = Constraint<EffectiveScalarMultivariateFunction,EffectiveNumber>; //!< <p/>
using ValidatedConstraint = Constraint<ValidatedScalarMultivariateFunction,ValidatedNumber>; //!< <p/>
using ApproximateConstraint = Constraint<ApproximateScalarMultivariateFunction,ApproximateNumber>; //!< <p/>
using ValidatedExactConstraint = Constraint<ValidatedScalarMultivariateFunction,ExactNumber>; //!< <p/>
using EffectiveExactConstraint = Constraint<EffectiveScalarMultivariateFunction,ExactNumber>; //!< <p/>

using EffectiveLowerConstraint = LowerConstraint<EffectiveScalarMultivariateFunction,EffectiveNumber>;
using ValidatedExactLowerConstraint = LowerConstraint<ValidatedScalarMultivariateFunction,ExactNumber>;
using ValidatedLowerConstraint = LowerConstraint<ValidatedScalarMultivariateFunction,ValidatedNumber>;

using EffectiveUpperConstraint = UpperConstraint<EffectiveScalarMultivariateFunction,EffectiveNumber>;
using ValidatedExactUpperConstraint = UpperConstraint<ValidatedScalarMultivariateFunction,ExactNumber>;
using ValidatedUpperConstraint = UpperConstraint<ValidatedScalarMultivariateFunction,ValidatedNumber>;

using EffectiveEqualityConstraint = EqualityConstraint<EffectiveScalarMultivariateFunction,EffectiveNumber>;
using ValidatedExactEqualityConstraint = EqualityConstraint<ValidatedScalarMultivariateFunction,ExactNumber>;
using ValidatedEqualityConstraint = EqualityConstraint<ValidatedScalarMultivariateFunction,ValidatedNumber>;
//!@}

inline EffectiveLowerConstraint operator<=(const EffectiveNumber& l, const EffectiveScalarMultivariateFunction& f) {
    return EffectiveLowerConstraint(l,f); }
inline EffectiveUpperConstraint operator>=(const EffectiveNumber& u, const EffectiveScalarMultivariateFunction& f) {
    return EffectiveUpperConstraint(f,u); }
inline EffectiveLowerConstraint operator>=(const EffectiveScalarMultivariateFunction& f, const EffectiveNumber& l) {
    return EffectiveLowerConstraint(l,f); }
inline EffectiveUpperConstraint operator<=(const EffectiveScalarMultivariateFunction& f, const EffectiveNumber& u) {
    return EffectiveUpperConstraint(f,u); }
inline EffectiveConstraint operator==(const EffectiveScalarMultivariateFunction& f, const EffectiveNumber& c) {
    return EffectiveConstraint(c,f,c); }

inline ValidatedExactLowerConstraint operator<=(const ExactNumber& l, const ValidatedScalarMultivariateFunction& f) {
    return ValidatedExactLowerConstraint(l,f); }
inline ValidatedExactUpperConstraint operator>=(const ExactNumber& u, const ValidatedScalarMultivariateFunction& f) {
    return ValidatedExactUpperConstraint(f,u); }
inline ValidatedExactLowerConstraint operator>=(const ValidatedScalarMultivariateFunction& f, const ExactNumber& l) {
    return ValidatedExactLowerConstraint(l,f); }
inline ValidatedExactUpperConstraint operator<=(const ValidatedScalarMultivariateFunction& f, const ExactNumber& u) {
    return ValidatedExactUpperConstraint(f,u); }
inline ValidatedExactConstraint operator==(const ValidatedScalarMultivariateFunction& f, const ExactNumber& c) {
    return ValidatedExactConstraint(c,f,c); }

inline ValidatedLowerConstraint operator<=(const ValidatedNumber& l, const ValidatedScalarMultivariateFunction& f) {
    return ValidatedLowerConstraint(l,f); }
inline ValidatedUpperConstraint operator>=(const ValidatedNumber& u, const ValidatedScalarMultivariateFunction& f) {
    return ValidatedUpperConstraint(f,u); }
inline ValidatedLowerConstraint operator>=(const ValidatedScalarMultivariateFunction& f, const ValidatedNumber& l) {
    return ValidatedLowerConstraint(l,f); }
inline ValidatedUpperConstraint operator<=(const ValidatedScalarMultivariateFunction& f, const ValidatedNumber& u) {
    return ValidatedUpperConstraint(f,u); }
inline ValidatedConstraint operator==(const ValidatedScalarMultivariateFunction& f, const ValidatedNumber& c) {
    return ValidatedConstraint(c,f,c); }

inline EffectiveConstraint operator<=(const EffectiveLowerConstraint& lc, const EffectiveNumber& u) {
    return EffectiveConstraint(lc.lower_bound(),lc.function(),u); }
inline ValidatedExactConstraint operator<=(const ValidatedExactLowerConstraint& lc, const ExactNumber& u) {
    return ValidatedExactConstraint(lc.lower_bound(),lc.function(),u); }
inline ValidatedConstraint operator<=(const ValidatedLowerConstraint& lc, const ValidatedNumber& u) {
    return ValidatedConstraint(lc.lower_bound(),lc.function(),u); }

inline EffectiveConstraint operator>=(const EffectiveUpperConstraint& uc, const EffectiveNumber& l) {
    return EffectiveConstraint(l,uc.function(),uc.upper_bound()); }
inline ValidatedExactConstraint operator>=(const ValidatedExactUpperConstraint& uc, const ExactNumber& l) {
    return ValidatedExactConstraint(l,uc.function(),uc.upper_bound()); }
inline ValidatedConstraint operator>=(const ValidatedUpperConstraint& uc, const ValidatedNumber& l) {
    return ValidatedConstraint(l,uc.function(),uc.upper_bound()); }

inline ValidatedExactConstraint operator<=(const ValidatedScalarMultivariateFunction& f1, const ValidatedScalarMultivariateFunction& f2) {
    return (f1-f2) <= ExactNumber(0); }
inline ValidatedExactConstraint operator>=(const ValidatedScalarMultivariateFunction& f1, const ValidatedScalarMultivariateFunction& f2) {
    return (f1-f2) >= ExactNumber(0); }

/*
template<class PF, class PR> inline LowerConstraint<ScalarMultivariateFunction<PF>,Number<PR>>
operator<=(const Number<PR>& l, ScalarMultivariateFunction<PF> const& f) {
    return LowerConstraint<ScalarMultivariateFunction<PF>,Number<PR>>(l,f); }
template<class PF, class PR> inline UpperConstraint<ScalarMultivariateFunction<PF>,Number<PR>>
operator<=(ScalarMultivariateFunction<PF> const& f, const Number<PR>& u) {
    return UpperConstraint<ScalarMultivariateFunction<PF>,Number<PR>>(f,u); }
template<class PF, class PR> inline UpperConstraint<ScalarMultivariateFunction<PF>,Number<PR>>
operator>=(const Number<PR>& u, ScalarMultivariateFunction<PF> const& f) {
    return UpperConstraint<ScalarMultivariateFunction<PF>,Number<PR>>(f,u); }
template<class PF, class PR> inline LowerConstraint<ScalarMultivariateFunction<PF>,Number<PR>>
operator>=(ScalarMultivariateFunction<PF> const& f, const Number<PR>& l) {
    return LowerConstraint<ScalarMultivariateFunction<PF>,Number<PR>>(l,f); }

template<class PF, class PR> inline Constraint<ScalarMultivariateFunction<PF>,Number<PR>>
operator<=(const LowerConstraint<ScalarMultivariateFunction<PF>,Number<PR>>& lc, Number<SelfType<PR>> const& u) {
    return Constraint<ScalarMultivariateFunction<PF>,Number<PR>>(lc.lower_bound(),lc.function(),u); }
template<class PF, class PR> inline Constraint<ScalarMultivariateFunction<PF>,Number<PR>>
operator<=(Number<SelfType<PR>> const& l, const UpperConstraint<ScalarMultivariateFunction<PF>,Number<PR>>& uc) {
    return Constraint<ScalarMultivariateFunction<PF>,Number<PR>>(l,uc.function(),uc.upper_bound()); }
template<class PF, class PR> inline Constraint<ScalarMultivariateFunction<PF>,Number<PR>>
operator>=(const UpperConstraint<ScalarMultivariateFunction<PF>,Number<PR>>& uc, Number<SelfType<PR>> const& l) {
    return Constraint<ScalarMultivariateFunction<PF>,Number<PR>>(l,uc.function(),uc.upper_bound()); }
template<class PF, class PR> inline Constraint<ScalarMultivariateFunction<PF>,Number<PR>>
operator>=(Number<SelfType<PR>> const& u, const LowerConstraint<ScalarMultivariateFunction<PF>,Number<PR>>& lc) {
    return Constraint<ScalarMultivariateFunction<PF>,Number<PR>>(lc.lower_bound(),lc.function(),u); }
*/


/*
inline EffectiveConstraint operator<=(const EffectiveNumber& l, const EffectiveScalarMultivariateFunction& f) {
    return EffectiveConstraint(l,f,+infty); }
inline EffectiveConstraint operator>=(const EffectiveNumber& u, const EffectiveScalarMultivariateFunction& f) {
    return EffectiveConstraint(-infty,f,u); }
inline EffectiveConstraint operator>=(const EffectiveScalarMultivariateFunction& f, const EffectiveNumber& l) {
    return EffectiveConstraint(l,f,+infty); }
inline EffectiveConstraint operator<=(const EffectiveScalarMultivariateFunction& f, const EffectiveNumber& u) {
    return EffectiveConstraint(-infty,f,u); }
inline EffectiveConstraint operator==(const EffectiveScalarMultivariateFunction& f, const EffectiveNumber& c) {
    return EffectiveConstraint(c,f,c); }

inline ValidatedExactConstraint operator<=(const ExactNumber& l, const ValidatedScalarMultivariateFunction& f) {
    return ValidatedExactConstraint(l,f,+infty); }
inline ValidatedExactConstraint operator>=(const ExactNumber& u, const ValidatedScalarMultivariateFunction& f) {
    return ValidatedExactConstraint(-infty,f,u); }
inline ValidatedExactConstraint operator>=(const ValidatedScalarMultivariateFunction& f, const ExactNumber& l) {
    return ValidatedExactConstraint(l,f,+infty); }
inline ValidatedExactConstraint operator<=(const ValidatedScalarMultivariateFunction& f, const ExactNumber& u) {
    return ValidatedExactConstraint(-infty,f,u); }
inline ValidatedExactConstraint operator==(const ValidatedScalarMultivariateFunction& f, const ExactNumber& c) {
    return ValidatedExactConstraint(c,f,c); }

inline ValidatedConstraint operator<=(const ValidatedNumber& l, const ValidatedScalarMultivariateFunction& f) {
    return ValidatedConstraint(l,f,+infty); }
inline ValidatedConstraint operator>=(const ValidatedNumber& u, const ValidatedScalarMultivariateFunction& f) {
    return ValidatedConstraint(-infty,f,u); }
inline ValidatedConstraint operator>=(const ValidatedScalarMultivariateFunction& f, const ValidatedNumber& l) {
    return ValidatedConstraint(l,f,+infty); }
inline ValidatedConstraint operator<=(const ValidatedScalarMultivariateFunction& f, const ValidatedNumber& u) {
    return ValidatedConstraint(-infty,f,u); }
inline ValidatedConstraint operator==(const ValidatedScalarMultivariateFunction& f, const ValidatedNumber& c) {
    return ValidatedConstraint(c,f,c); }


inline EffectiveExactConstraint operator<=(const EffectiveExactConstraint& nc, const ExactNumber& c) {
    ARIADNE_ASSERT(decide(nc.upper_bound()==infty));
    return EffectiveExactConstraint(nc.lower_bound(),nc.function(),c);
}
inline EffectiveConstraint operator<=(const EffectiveConstraint& nc, const EffectiveNumber& c) {
    ARIADNE_ASSERT(decide(nc.upper_bound()==infty));
    return EffectiveConstraint(nc.lower_bound(),nc.function(),c);
}
inline ValidatedExactConstraint operator<=(const ValidatedExactConstraint& nc, const ExactNumber& c) {
    ARIADNE_ASSERT(decide(nc.upper_bound()==infty));
    return ValidatedExactConstraint(nc.lower_bound(),nc.function(),c);
}
inline ValidatedConstraint operator<=(const ValidatedConstraint& nc, const ValidatedNumber& c) {
    ARIADNE_ASSERT(decide(nc.upper_bound()==infty));
    return ValidatedConstraint(nc.lower_bound(),nc.function(),c);
}


inline ValidatedExactConstraint operator<=(const ValidatedScalarMultivariateFunction& f1, const ValidatedScalarMultivariateFunction& f2) {
    return (f1-f2) <= ExactNumber(0); }
inline ValidatedExactConstraint operator>=(const ValidatedScalarMultivariateFunction& f1, const ValidatedScalarMultivariateFunction& f2) {
    return (f1-f2) >= ExactNumber(0); }
*/


/*
inline EffectiveConstraint operator<=(const ExactNumber& l, const EffectiveScalarMultivariateFunction& f) {
    return Constraint<ScalarMultivariateFunction<PF>,ExactNumber>(l,f,+infty); }
inline EffectiveConstraint operator>=(const EffectiveNumber& u, const EffectiveScalarMultivariateFunction& f) {
    return EffectiveConstraint(-infty,f,u); }
inline EffectiveConstraint operator>=(const EffectiveScalarMultivariateFunction& f, const EffectiveNumber& l) {
    return EffectiveConstraint(l,f,+infty); }
inline EffectiveConstraint operator<=(const EffectiveScalarMultivariateFunction& f, const EffectiveNumber& u) {
    return EffectiveConstraint(-infty,f,u); }


inline EffectiveConstraint operator<=(const EffectiveNumber& l, const EffectiveScalarMultivariateFunction& f) {
    return EffectiveConstraint(l,f,+infty); }
inline EffectiveConstraint operator>=(const EffectiveNumber& u, const EffectiveScalarMultivariateFunction& f) {
    return EffectiveConstraint(-infty,f,u); }
inline EffectiveConstraint operator>=(const EffectiveScalarMultivariateFunction& f, const EffectiveNumber& l) {
    return EffectiveConstraint(l,f,+infty); }
inline EffectiveConstraint operator<=(const EffectiveScalarMultivariateFunction& f, const EffectiveNumber& u) {
    return EffectiveConstraint(-infty,f,u); }

inline EffectiveExactConstraint operator<=(const ExactNumber& l, const EffectiveScalarMultivariateFunction& f) {
    return EffectiveExactConstraint(l,f,+infty); }
inline EffectiveExactConstraint operator>=(const ExactNumber& u, const EffectiveScalarMultivariateFunction& f) {
    return EffectiveExactConstraint(-infty,f,u); }
inline EffectiveExactConstraint operator>=(const EffectiveScalarMultivariateFunction& f, const ExactNumber& l) {
    return EffectiveExactConstraint(l,f,+infty); }
inline EffectiveExactConstraint operator<=(const EffectiveScalarMultivariateFunction& f, const ExactNumber& u) {
    return EffectiveExactConstraint(-infty,f,u); }

inline ValidatedExactConstraint operator<=(const ExactNumber& l, const ValidatedScalarMultivariateFunction& f) {
    return ValidatedExactConstraint(l,f,+infty); }
inline ValidatedExactConstraint operator>=(const ExactNumber& u, const ValidatedScalarMultivariateFunction& f) {
    return ValidatedExactConstraint(-infty,f,u); }
inline ValidatedExactConstraint operator>=(const ValidatedScalarMultivariateFunction& f, const ExactNumber& l) {
    return ValidatedExactConstraint(l,f,+infty); }
inline ValidatedExactConstraint operator<=(const ValidatedScalarMultivariateFunction& f, const ExactNumber& u) {
    return ValidatedExactConstraint(-infty,f,u); }
*/

/*

template<class PF, class PR> inline Constraint<ScalarMultivariateFunction<PF>,Number<PR>>
operator<=(const Number<PR>& l, ScalarMultivariateFunction<PF> const& f) {
    return Constraint<ScalarMultivariateFunction<PF>,Number<PR>>(l,f,infty); }

template<class PF, class PR> inline Constraint<ScalarMultivariateFunction<PF>,Number<PR>>
operator<=(const Constraint<ScalarMultivariateFunction<PF>,Number<PR>>& lc, Number<SelfType<PR>> const& u) {
    return Constraint<ScalarMultivariateFunction<PF>,Number<PR>>(lc.lower_bound(),lc.function(),min(lc.upper_bound(),u)); }

template<class PF, class PR> inline Constraint<ScalarMultivariateFunction<PF>,Number<PR>>
operator==(const ScalarMultivariateFunction<SelfType<PF>>& f, const Number<PR>& c) {
    return Constraint<ScalarMultivariateFunction<PF>,Number<PR>>(f,c);
}

template<class PF, class PR> inline Constraint<ScalarMultivariateFunction<PF>,Number<PR>>
operator<=(ScalarMultivariateFunction<PF> const& f, const Number<PR>& u) {
    return Constraint<ScalarMultivariateFunction<PF>,Number<PR>>(-infty,f,u); }

template<class PF, class PR> inline Constraint<ScalarMultivariateFunction<PF>,Number<PR>>
operator>=(const Number<PR>& u, ScalarMultivariateFunction<PF> const& f) {
    return Constraint<ScalarMultivariateFunction<PF>,Number<PR>>(-infty,f,u); }

template<class PF, class PR> inline Constraint<ScalarMultivariateFunction<PF>,Number<PR>>
operator>=(ScalarMultivariateFunction<PF> const& f, const Number<PR>& l) {
    return Constraint<ScalarMultivariateFunction<PF>,Number<PR>>(l,f,infty); }

*/

/*

inline EffectiveConstraint operator<=(const EffectiveNumber& c, const EffectiveScalarMultivariateFunction& f) {
    return EffectiveConstraint(c,f,infinity);
}

inline EffectiveConstraint operator>=(const EffectiveNumber& c, const EffectiveScalarMultivariateFunction& f) {
    return EffectiveConstraint(-infinity,f,c);
}

inline EffectiveConstraint operator<=(const EffectiveScalarMultivariateFunction& f, const EffectiveNumber& c) {
    return EffectiveConstraint(-infinity,f,c);
}

inline EffectiveConstraint operator>=(const EffectiveScalarMultivariateFunction& f, const EffectiveNumber& c) {
    return EffectiveConstraint(c,f,infinity);
}

inline EffectiveConstraint operator==(const EffectiveScalarMultivariateFunction& f, const EffectiveNumber& c) {
    return EffectiveConstraint(f,c);
}

inline EffectiveConstraint operator<=(const EffectiveScalarMultivariateFunction& f, double c) {
    return f <= Dyadic(c);
}

inline EffectiveConstraint operator>=(const EffectiveScalarMultivariateFunction& f, double c) {
    return f >= Dyadic(c);
}

inline EffectiveConstraint operator==(const EffectiveScalarMultivariateFunction& f, double c) {
    return f == Dyadic(c);
}

inline EffectiveConstraint operator<=(const EffectiveConstraint& nc, const EffectiveNumber& c) {
    ARIADNE_ASSERT(decide(nc.upper_bound()==infty));
    return EffectiveConstraint(nc.lower_bound(),nc.function(),c);
}


inline ValidatedExactConstraint operator<=(const ExactNumber& c, const ValidatedScalarMultivariateFunction& f) {
    return ValidatedExactConstraint(c,f,ExactNumber(+infty));
}

inline ValidatedExactConstraint operator<=(const ValidatedScalarMultivariateFunction& f, const ExactNumber& c) {
    return ValidatedExactConstraint(-infty,f,c);
}

inline ValidatedExactConstraint operator>=(const ValidatedScalarMultivariateFunction& f, const ExactNumber& c) {
    return ValidatedExactConstraint(c,f,+infty);
}

inline ValidatedExactConstraint operator==(const ValidatedScalarMultivariateFunction& f, const ExactNumber& c) {
    return ValidatedExactConstraint(c,f,c);
}

inline ValidatedExactConstraint operator<=(const ValidatedExactConstraint& nc, const ExactNumber& c) {
    ARIADNE_ASSERT(nc.upper_bound()==infty);
    return ValidatedExactConstraint(nc.lower_bound(),nc.function(),c);
}


inline ValidatedExactConstraint operator<=(const ValidatedScalarMultivariateFunction& f1, const ValidatedScalarMultivariateFunction& f2) {
    return (f1-f2) <= ExactNumber(0);
}

inline ValidatedExactConstraint operator>=(const ValidatedScalarMultivariateFunction& f1, const ValidatedScalarMultivariateFunction& f2) {
    return (f1-f2) >= ExactNumber(0);
}


inline ValidatedConstraint operator<=(const ValidatedNumber& c, const ValidatedScalarMultivariateFunction& f) {
    return ValidatedConstraint(c,f,ExactNumber(+infty));
}

inline ValidatedConstraint operator<=(const ValidatedScalarMultivariateFunction& f, const ValidatedNumber& c) {
    return ValidatedConstraint(ExactNumber(-infty),f,c);
}

inline ValidatedConstraint operator<=(const ValidatedConstraint& nc, const ValidatedNumber& c) {
    ARIADNE_ASSERT(decide(nc.upper_bound()==infty));
    return ValidatedConstraint(nc.lower_bound(),nc.function(),c);
}

*/

} //namespace Ariadne

#endif
