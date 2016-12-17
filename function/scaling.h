/***************************************************************************
 *            scaling.h
 *
 *  Copyright 2008-15  Pieter Collins
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

/*! \file scaling.h
 *  \brief Scaling functions.
 */

#ifndef ARIADNE_SCALING_H
#define ARIADNE_SCALING_H

#include "numeric/dyadic.h"
#include "geometry/interval.h"
#include "geometry/box.h"

namespace Ariadne {

inline ApproximateNumericType med_apprx(ExactIntervalType const& ivl) {
    return ApproximateNumericType(hlf_exact(add_approx(ivl.lower().raw(),ivl.upper().raw())));
}

inline ApproximateNumericType rad_apprx(ExactIntervalType const& ivl) {
    return ApproximateNumericType(hlf_exact(sub_approx(ivl.upper().raw(),ivl.lower().raw())));
}

inline ValidatedNumericType med_val(ExactIntervalType const& ivl) {
    return hlf(ivl.lower()+ivl.upper());
}

inline ValidatedNumericType rad_val(ExactIntervalType const& ivl) {
    return hlf(ivl.upper()-ivl.lower());
}

inline Dyadic med(ExactIntervalType const& ivl) {
    return hlf(add( Dyadic(ivl.lower().raw()), Dyadic(ivl.upper().raw()) ));
}

inline Dyadic rad(ExactIntervalType const& ivl) {
    return hlf(sub( Dyadic(ivl.upper().raw()), Dyadic(ivl.lower().raw()) ));
}

template<class X> class Formula;
template<class X> class Algebra;
template<class X> Formula<X> unscale(Formula<X> x, const ExactIntervalType& d) { assert(false); }
template<class X> Algebra<X> unscale(Algebra<X> x, const ExactIntervalType& d) { assert(false); }

template<class T, EnableIf<IsSame<Paradigm<T>,ApproximateTag>> =dummy>
inline T unscale(T x, const ExactIntervalType& d) {
    return (std::move(x)-med(d))/rad(d);
}

template<class T, EnableIf<IsStronger<Paradigm<T>,ValidatedTag>> =dummy>
inline T unscale(T x, const ExactIntervalType& d) {
    if(d.lower()==d.upper()) { return std::move(x)*Dyadic(0); }
    return (std::move(x)-med(d))/rad(d);
}

template<class X> Vector<X> unscale(const Vector<X>& x, const ExactBoxType& d) {
    Vector<X> r(x);
    for(SizeType i=0; i!=r.size(); ++i) {
        r[i]=unscale(x[i],d[i]);
    }
    return r;
}

class Scaling {
    ExactIntervalType _codom;
  public:
    Scaling(ExactIntervalType codom) : _codom(codom) { }
    UnitIntervalType domain() const { return UnitIntervalType(); }
    ExactIntervalType codomain() const { return _codom; }
    template<class X> X operator() (X) const;
};

class VectorScaling {
    Box<ExactIntervalType> _codom;
  public:
    VectorScaling(Box<ExactIntervalType> codom) : _codom(codom) { }
    SizeType size() const { return _codom.dimension(); }
    Scaling operator[] (SizeType i) const { return Scaling(_codom[i]); }
    Box<ExactIntervalType> const& codomain() const { return _codom; }
    template<class X> Vector<X> operator() (Vector<X> const&) const;
};

class Unscaling {
    ExactIntervalType _dom;
  public:
    Unscaling(ExactIntervalType dom) : _dom(dom) { }
    ExactIntervalType domain() const { return _dom; }
    UnitIntervalType codomain() const { return UnitIntervalType(); }
    template<class X> X operator() (X) const;
};

class VectorUnscaling {
    Box<ExactIntervalType> _dom;
  public:
    VectorUnscaling(Box<ExactIntervalType> dom) : _dom(dom) { }
    SizeType size() const { return _dom.dimension(); }
    Unscaling operator[] (SizeType i) const { return Unscaling(_dom[i]); }
    Box<ExactIntervalType> const& domain() const { return _dom; }
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

#endif // ARIADNE_SCALING_H
