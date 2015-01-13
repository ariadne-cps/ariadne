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


class Scaling {
    ExactInterval _codom;
  public:
    Scaling(ExactInterval codom) : _codom(codom) { }
    UnitInterval domain() const { return UnitInterval(); }
    ExactInterval codomain() const { return _codom; }
    template<class X> X operator() (X const&) const;
};

class VectorScaling {
    Box<ExactInterval> _codom;
  public:
    VectorScaling(Box<ExactInterval> codom) : _codom(codom) { }
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
    template<class X> X operator() (X const&) const;
};

class VectorUnscaling {
    Box<ExactInterval> _dom;
  public:
    VectorUnscaling(Box<ExactInterval> dom) : _dom(dom) { }
    Unscaling operator[] (SizeType i) const { return Unscaling(_dom[i]); }
    Box<ExactInterval> const& domain() const { return _dom; }
    template<class X> Vector<X> operator() (Vector<X> const&) const;
};

} // namespace Ariadne

#endif // ARIADNE_SCALING_H
