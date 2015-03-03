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

#include "geometry/interval.h"
#include "geometry/box.h"

namespace Ariadne {

inline ApproximateNumericType med_apprx(ExactInterval const& ivl) {
    return ApproximateNumericType(half_exact(add_approx(ivl.lower().raw(),ivl.upper().raw())));
}

inline ApproximateNumericType rad_apprx(ExactInterval const& ivl) {
    return ApproximateNumericType(half_exact(sub_approx(ivl.upper().raw(),ivl.lower().raw())));
}

inline ValidatedNumericType med_val(ExactInterval const& ivl) {
    return half(ivl.lower()+ivl.upper());
}

inline ValidatedNumericType rad_val(ExactInterval const& ivl) {
    return half(ivl.upper()-ivl.lower());
}


template<class T, EnableIf<IsSame<Paradigm<T>,Approximate>> =dummy>
inline T unscale(T x, const ExactInterval& d) {
    ApproximateNumericType c(med_apprx(d));
    ApproximateNumericType r(rad_apprx(d));
    return (std::move(x)-c)/r;
}

template<class T, EnableIf<IsStronger<Paradigm<T>,Validated>> =dummy>
inline T unscale(T x, const ExactInterval& d) {
    ValidatedNumericType c(med_val(d));
    if(d.lower()==d.upper()) { c=0; return std::move(x)*c; }
    ValidatedNumericType r(rad_val(d));
    return (std::move(x)-c)/r;
}

template<class X> Vector<X> unscale(const Vector<X>& x, const ExactBox& d) {
    Vector<X> r(x);
    for(SizeType i=0; i!=r.size(); ++i) {
        r[i]=unscale(x[i],d[i]);
    }
    return r;
}

class Scaling {
    ExactInterval _codom;
  public:
    Scaling(ExactInterval codom) : _codom(codom) { }
    UnitInterval domain() const { return UnitInterval(); }
    ExactInterval codomain() const { return _codom; }
    template<class X> X operator() (X) const;
};

class VectorScaling {
    Box<ExactInterval> _codom;
  public:
    VectorScaling(Box<ExactInterval> codom) : _codom(codom) { }
    SizeType size() const { return _codom.dimension(); }
    Scaling operator[] (SizeType i) const { return Scaling(_codom[i]); }
    Box<ExactInterval> const& codomain() const { return _codom; }
    template<class X> Vector<X> operator() (Vector<X> const&) const;
};

class Unscaling {
    ExactInterval _dom;
  public:
    Unscaling(ExactInterval dom) : _dom(dom) { }
    ExactInterval domain() const { return _dom; }
    UnitInterval codomain() const { return UnitInterval(); }
    template<class X> X operator() (X) const;
};

class VectorUnscaling {
    Box<ExactInterval> _dom;
  public:
    VectorUnscaling(Box<ExactInterval> dom) : _dom(dom) { }
    SizeType size() const { return _dom.dimension(); }
    Unscaling operator[] (SizeType i) const { return Unscaling(_dom[i]); }
    Box<ExactInterval> const& domain() const { return _dom; }
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
