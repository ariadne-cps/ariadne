/***************************************************************************
 *            attribute.hpp
 *
 *  Copyright 2011-17  Pieter Collins
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

/*! \file attribute.hpp
 *  \brief Generic attributes for named parameters.
 */

#ifndef ARIADNE_ATTRIBUTE_HPP
#define ARIADNE_ATTRIBUTE_HPP

namespace Ariadne {

template<class T> class Generator {
  public:
    inline T operator=(const typename T::Type& v) const;
};

template<class V> class Attribute {
    template<class T> friend class Generator;
  private:
    V _v;
  protected:
//    Attribute(const V& v) : _v(v) { }
  public:
//    explicit Attribute(const V& v) : _v(v) { }
    Attribute(const V& v) : _v(v) { }
  public:
    typedef V Type;
    operator V() const { return this->_v; }
};

template<class T> inline T Generator<T>::operator=(const typename T::Type& v) const {
    Attribute<typename T::Type> attr(v); return static_cast<const T&>(attr); }

struct MaximumError : Attribute<double> { using Attribute<double>::Attribute; };
struct SweepThreshold : Attribute<double> { using Attribute<double>::Attribute; };
struct MaximumNumericTypeOfSteps : Attribute<uint> { using Attribute<uint>::Attribute; };

static const Generator<MaximumError> maximum_error = Generator<MaximumError>();
static const Generator<SweepThreshold> sweep_threshold = Generator<SweepThreshold>();
static const Generator<MaximumNumericTypeOfSteps> maximum_number_of_steps = Generator<MaximumNumericTypeOfSteps>();

struct Capacity : Attribute<SizeType> { };
static const Generator<Capacity> capacity = Generator<Capacity>();
struct Size : Attribute<SizeType> { };
static const Generator<Size> size = Generator<Size>();
struct ResultSize : Attribute<SizeType> { };
static const Generator<ResultSize> result_size = Generator<ResultSize>();
struct ArgumentSize : Attribute<SizeType> { };
static const Generator<ArgumentSize> argument_size = Generator<ArgumentSize>();
struct Degree : Attribute<DegreeType> { };
static const Generator<Degree> degree = Generator<Degree>();


} // namespace Ariadne

#endif /* ARIADNE_ATTRIBUTE_HPP */
