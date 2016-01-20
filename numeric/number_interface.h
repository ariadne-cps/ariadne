/***************************************************************************
 *            numeric/number_interface.h
 *
 *  Copyright 2013-14  Pieter Collins
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

/*! \file numeric/number_interface.h
 *  \brief
 */



#ifndef ARIADNE_NUMBER_INTERFACE
#define ARIADNE_NUMBER_INTERFACE

#include "utility/clonable.h"
#include "utility/writable.h"

#include "number.decl.h"
#include "float.decl.h"

namespace Ariadne {

/************ Number *********************************************************/

class NumberInterface;
class UnaryOperatorInterface;



class NumberInterface
    : public virtual WritableInterface
    , public virtual ClonableInterface
{
    template<class N> friend class NumberWrapper;
    friend class Handle<NumberInterface>;
  public:
    virtual ~NumberInterface() { }
  public:
    virtual NumberInterface* _copy() const = 0;
    virtual NumberInterface* _move() = 0;
    virtual NumberInterface* _pos() const = 0;
    virtual NumberInterface* _neg() const = 0;
    virtual NumberInterface* _sqr() const = 0;
    virtual NumberInterface* _rec() const = 0;
    virtual NumberInterface* _pow(Int n) const = 0;
    virtual NumberInterface* _add(NumberInterface const& y) const = 0;
    virtual NumberInterface* _sub(NumberInterface const& y) const = 0;
    virtual NumberInterface* _mul(NumberInterface const& y) const = 0;
    virtual NumberInterface* _div(NumberInterface const& y) const = 0;
    virtual NumberInterface* _sqrt() const = 0;
    virtual NumberInterface* _exp() const = 0;
    virtual NumberInterface* _log() const = 0;
    virtual NumberInterface* _sin() const = 0;
    virtual NumberInterface* _cos() const = 0;
    virtual NumberInterface* _tan() const = 0;
    virtual NumberInterface* _atan() const = 0;
    virtual NumberInterface* _abs() const = 0;
    virtual NumberInterface* _max(NumberInterface const& y) const = 0;
    virtual NumberInterface* _min(NumberInterface const& y) const = 0;
    virtual NumberInterface* _apply(UnaryOperatorInterface const& o) const = 0;

    virtual Float64Ball _get(MetricTag, Precision64) const = 0;
    virtual Float64Bounds _get(BoundedTag, Precision64) const = 0;
    virtual Float64UpperBound _get(UpperTag, Precision64) const = 0;
    virtual Float64LowerBound _get(LowerTag, Precision64) const = 0;
    virtual Float64Approximation _get(ApproximateTag, Precision64) const = 0;

    virtual FloatMPBall _get(MetricTag, PrecisionMP) const = 0;
    virtual FloatMPBounds _get(BoundedTag, PrecisionMP) const = 0;
    virtual FloatMPUpperBound _get(UpperTag, PrecisionMP) const = 0;
    virtual FloatMPLowerBound _get(LowerTag, PrecisionMP) const = 0;
    virtual FloatMPApproximation _get(ApproximateTag, PrecisionMP) const = 0;

    virtual ParadigmCode _paradigm() const = 0;
    virtual String _class_name() const = 0;

    virtual LogicalValue _equals(NumberInterface const& y) const = 0;
    virtual LogicalValue _less(NumberInterface const& y) const = 0;
};


} // namespace Ariadne

#endif /* ARIADNE_NUMBER_INTERFACE */
