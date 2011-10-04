/***************************************************************************
 *            attribute.h
 *
 *  Copyright  2011  Pieter Collins
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

/*! \file attribute.h
 *  \brief Generic attributes for named parameters.
 */

#ifndef ARIADNE_ATTRIBUTE_H
#define ARIADNE_ATTRIBUTE_H

namespace Ariadne {

template<class V> class Attribute {
    V _v;
  protected:
    Attribute(const V& v) : _v(v) { }
  public:
    typedef V Type;
    operator V() const { return this->_v; }
};

template<class T> class Generator {
  public:
    T operator=(const typename T::Type& v) const { return T(v); }
};

struct MaximumError : Attribute<double> { MaximumError(double v) : Attribute(v) { } };
struct SweepThreshold : Attribute<double> { SweepThreshold(double v) : Attribute(v) { } };
struct MaximumNumberOfSteps : Attribute<double> { MaximumNumberOfSteps(double v) : Attribute(v) { } };

static const Generator<MaximumError> maximum_error = Generator<MaximumError>();
static const Generator<SweepThreshold> sweep_threshold = Generator<SweepThreshold>();
static const Generator<MaximumNumberOfSteps> maximum_number_of_steps = Generator<MaximumNumberOfSteps>();

} // namespace Ariadne

#endif /* ARIADNE_ATTRIBUTE_H */
