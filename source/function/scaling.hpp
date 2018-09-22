/***************************************************************************
 *            scaling.hpp
 *
 *  Copyright 2008-17  Pieter Collins
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

/*! \file scaling.hpp
 *  \brief Scaling functions.
 */

#ifndef ARIADNE_SCALING_HPP
#define ARIADNE_SCALING_HPP

#include "../numeric/dyadic.hpp"
#include "../function/domain.hpp"

namespace Ariadne {

struct UnscalingException : std::runtime_error {
    IntervalDomainType domain;
    inline UnscalingException(String msg, IntervalDomainType dom) : std::runtime_error(msg+to_str(dom)), domain(dom) { }
};

inline ApproximateNumericType med_apprx(IntervalDomainType const& ivl) {
    return ApproximateNumericType(hlf(add(approx,ivl.lower().raw(),ivl.upper().raw())));
}

inline ApproximateNumericType rad_apprx(IntervalDomainType const& ivl) {
    return ApproximateNumericType(hlf(sub(approx,ivl.upper().raw(),ivl.lower().raw())));
}

inline ValidatedNumericType med_val(IntervalDomainType const& ivl) {
    return hlf(ivl.lower()+ivl.upper());
}

inline ValidatedNumericType rad_val(IntervalDomainType const& ivl) {
    return hlf(ivl.upper()-ivl.lower());
}

inline Dyadic med(IntervalDomainType const& ivl) {
    return hlf(add( Dyadic(ivl.lower().raw()), Dyadic(ivl.upper().raw()) ));
}

inline Dyadic rad(IntervalDomainType const& ivl) {
    return hlf(sub( Dyadic(ivl.upper().raw()), Dyadic(ivl.lower().raw()) ));
}

template<class T> inline
T unscale(T x, const IntervalDomainType& d) {
    if(d.lower()==d.upper()) {
        return std::move(x)*0;
        // throw UnscalingException{"Cannot unscale "+to_str(x)+" over empty interval ",d};
    } else {
        return (std::move(x)-med(d))/rad(d);
    }
}

template<class X> Vector<X> unscale(const Vector<X>& x, const BoxDomainType& d) {
    Vector<X> r(x);
    for(SizeType i=0; i!=r.size(); ++i) {
        r[i]=unscale(x[i],d[i]);
    }
    return r;
}

class Scaling {
    IntervalDomainType _codom;
  public:
    Scaling(IntervalDomainType codom) : _codom(codom) { }
    UnitInterval domain() const { return UnitInterval(); }
    IntervalDomainType codomain() const { return _codom; }
    template<class X> X operator() (X) const;
};

class VectorScaling {
    Box<IntervalDomainType> _codom;
  public:
    VectorScaling(Box<IntervalDomainType> codom) : _codom(codom) { }
    SizeType size() const { return _codom.dimension(); }
    Scaling operator[] (SizeType i) const { return Scaling(_codom[i]); }
    Box<IntervalDomainType> const& codomain() const { return _codom; }
    template<class X> Vector<X> operator() (Vector<X> const&) const;
};

class Unscaling {
    IntervalDomainType _dom;
  public:
    Unscaling(IntervalDomainType dom) : _dom(dom) { }
    IntervalDomainType domain() const { return _dom; }
    UnitInterval codomain() const { return UnitInterval(); }
    template<class X> X operator() (X) const;
};

class VectorUnscaling {
    Box<IntervalDomainType> _dom;
  public:
    VectorUnscaling(Box<IntervalDomainType> dom) : _dom(dom) { }
    SizeType size() const { return _dom.dimension(); }
    Unscaling operator[] (SizeType i) const { return Unscaling(_dom[i]); }
    Box<IntervalDomainType> const& domain() const { return _dom; }
    template<class X> Vector<X> operator() (Vector<X> const&) const;
};

template<class X> X Scaling::operator() (X x) const {
    auto r=_codom.radius(); auto c=_codom.midpoint();
    return x*r+c;
}

template<class X> X Unscaling::operator() (X x) const {
    return unscale(std::move(x),this->_dom);
    auto r=_dom.radius(); auto c=_dom.midpoint();
    return (x-c)/r;
}

} // namespace Ariadne

#endif // ARIADNE_SCALING_HPP
