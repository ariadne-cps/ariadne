/***************************************************************************
 *            domain.h
 *
 *  Copyright 2008-13  Alberto Casagrande, Pieter Collins
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

/*! \file domain.h
 *  \brief Interval and box domains for functions.
 */

#ifndef ARIADNE_DOMAIN_H
#define ARIADNE_DOMAIN_H

#include "geometry/interval.h"
#include "geometry/box.h"

namespace Ariadne {

class RealDomain : public IntervalDomainType {
  public:
    RealDomain() : IntervalDomainType(-inf,+inf) { }
};

class EuclideanDomain : public BoxDomainType {
  public:
    EuclideanDomain(SizeType n) : BoxDomainType(n,RealDomain()) { }
};

using IntervalDomainType = ExactIntervalType;
using BoxDomainType = ExactBoxType;

inline SizeOne dimension(IntervalDomainType dom) { return SizeOne(); }
inline SizeType dimension(BoxDomainType dom) { return dom.dimension(); }


class UnitInterval;
typedef Box<UnitInterval> UnitBox;
enum class SplitPart : char;

template<class S> struct ElementTraits;
template<class S, class X> using ElementType = typename ElementTraits<S>::template Type<X>;
template<> struct ElementTraits<IntervalDomainType> { template<class X> using Type=Scalar<X>; };
template<> struct ElementTraits<BoxDomainType> { template<class X> using Type=Vector<X>; };

} // namespace Ariadne


#endif /* ARIADNE_DOMAIN_H */
