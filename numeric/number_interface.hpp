/***************************************************************************
 *            numeric/number_interface.h
 *
 *  Copyright 2013-17  Pieter Collins
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
#include "operators.h"

namespace Ariadne {

/************ Number *********************************************************/

class NumberInterface;
class UnaryOperatorInterface;
template<class... YS> struct Aware;
template<class X> class NumberWrapper;



class NumberInterface
    : public std::enable_shared_from_this<NumberInterface>
    , public virtual WritableInterface
    , public virtual ClonableInterface
{
    template<class X> friend class NumberWrapper;
    friend class Handle<NumberInterface>;
  public:
    virtual ~NumberInterface() { }
  public:
    virtual NumberInterface* _copy() const = 0;
    virtual NumberInterface* _move() = 0;

    virtual NumberInterface* _apply(Add op, NumberInterface const* y) const = 0;
    virtual NumberInterface* _apply(Sub op, NumberInterface const* y) const = 0;
    virtual NumberInterface* _apply(Mul op, NumberInterface const* y) const = 0;
    virtual NumberInterface* _apply(Div op, NumberInterface const* y) const = 0;
    virtual NumberInterface* _rapply(Add op, NumberInterface const* y) const = 0;
    virtual NumberInterface* _rapply(Sub op, NumberInterface const* y) const = 0;
    virtual NumberInterface* _rapply(Mul op, NumberInterface const* y) const = 0;
    virtual NumberInterface* _rapply(Div op, NumberInterface const* y) const = 0;

    virtual NumberInterface* _apply(Pos op) const = 0;
    virtual NumberInterface* _apply(Neg op) const = 0;
    virtual NumberInterface* _apply(Sqr op) const = 0;
    virtual NumberInterface* _apply(Rec op) const = 0;
    virtual NumberInterface* _apply(Pow op, Int n) const = 0;
    virtual NumberInterface* _apply(Sqrt op) const = 0;
    virtual NumberInterface* _apply(Exp op) const = 0;
    virtual NumberInterface* _apply(Log op) const = 0;
    virtual NumberInterface* _apply(Sin op) const = 0;
    virtual NumberInterface* _apply(Cos op) const = 0;
    virtual NumberInterface* _apply(Tan op) const = 0;
    virtual NumberInterface* _apply(Atan op) const = 0;

    virtual NumberInterface* _apply(Abs op) const = 0;
    virtual NumberInterface* _apply(Max op, NumberInterface const* y) const = 0;
    virtual NumberInterface* _apply(Min op, NumberInterface const* y) const = 0;
    virtual NumberInterface* _rapply(Max op, NumberInterface const* y) const = 0;
    virtual NumberInterface* _rapply(Min op, NumberInterface const* y) const = 0;

    virtual Rational _get_q() const = 0;

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
