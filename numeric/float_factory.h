/***************************************************************************
 *            float_factory.h
 *
 *  Copyright 2008-16  Pieter Collins
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

/*! \file float_factory.h
 *  \brief Factories for creating floating-point numbers.
 */

#ifndef ARIADNE_FLOAT_FACTORY_H
#define ARIADNE_FLOAT_FACTORY_H

#include "float-user.h"

namespace Ariadne {

template<class PR> class FloatFactory {
    PR _pr;
  public:
    typedef PR PrecisionType;
    FloatFactory(PR const& pr) : _pr(pr) { }
    PR precision() const { return this->_pr; }
  public:
    FloatApproximation<PR> create(Number<ApproximateTag> const& y) { return FloatApproximation<PR>(y,_pr); }
    FloatLowerBound<PR> create(Number<ValidatedLowerTag> const& y) { return FloatLowerBound<PR>(y,_pr); }
    FloatUpperBound<PR> create(Number<ValidatedUpperTag> const& y) { return FloatUpperBound<PR>(y,_pr); }
    FloatBounds<PR> create(Number<ValidatedTag> const& y) { return FloatBounds<PR>(y,_pr); }
    FloatBounds<PR> create(Number<EffectiveTag> const& y) { return FloatBounds<PR>(y,_pr); }
    FloatBounds<PR> create(Number<ExactTag> const& y) { return FloatBounds<PR>(y,_pr); }
    FloatBounds<PR> create(Real const& y) { return FloatBounds<PR>(y,_pr); }
    FloatBounds<PR> create(Rational const& y) { return FloatBounds<PR>(y,_pr); }
    FloatBounds<PR> create(Dyadic const& y) { return FloatBounds<PR>(y,_pr); }
    FloatBounds<PR> create(Integer const& y) { return FloatBounds<PR>(y,_pr); }
    FloatValue<PR> create(Dyadic const& y, ExactTag) { return FloatValue<PR>(y,_pr); }
    FloatValue<PR> create(Integer const& y, ExactTag) { return FloatValue<PR>(y,_pr); }
    template<class N, EnableIf<IsSignedIntegral<N>> =dummy> FloatValue<PR> create(N const& y) { return FloatValue<PR>(y,_pr); }
    template<class M, EnableIf<IsUnsignedIntegral<M>> =dummy> PositiveFloatValue<PR> create(M const& y) { return PositiveFloatValue<PR>(y,_pr); }
    template<class D, EnableIf<IsFloatingPoint<D>> =dummy> FloatApproximation<PR> create(D const& y) { return FloatApproximation<PR>(RawFloat<PR>(y,_pr)); }
};
template<class Y, class PR> using ConcreteType = decltype(declval<FloatFactory<PR>>().create(declval<Y>()));

template<class PR> inline FloatFactory<PR> float_factory(PR pr) { return FloatFactory<PR>(pr); }
template<class PR> inline FloatFactory<PR> factory(FloatApproximation<PR> const& flt) { return FloatFactory<PR>(flt.precision()); }
template<class PR> inline FloatFactory<PR> factory(FloatLowerBound<PR> const& flt) { return FloatFactory<PR>(flt.precision()); }
template<class PR> inline FloatFactory<PR> factory(FloatUpperBound<PR> const& flt) { return FloatFactory<PR>(flt.precision()); }
template<class PR> inline FloatFactory<PR> factory(FloatBounds<PR> const& flt) { return FloatFactory<PR>(flt.precision()); }
template<class PR> inline FloatFactory<PR> factory(FloatBall<PR> const& flt) { return FloatFactory<PR>(flt.precision()); }
template<class PR> inline FloatFactory<PR> factory(FloatValue<PR> const& flt) { return FloatFactory<PR>(flt.precision()); }

template<class Y, class PR> inline decltype(auto) make_float(Y const& y, PR pr) { return float_factory(pr).create(y); }

} // namespace Ariadne

#endif
