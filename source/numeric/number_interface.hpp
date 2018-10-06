/***************************************************************************
 *            numeric/number_interface.hpp
 *
 *  Copyright 2013-17  Pieter Collins
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

/*! \file numeric/number_interface.hpp
 *  \brief
 */



#ifndef ARIADNE_NUMBER_INTERFACE
#define ARIADNE_NUMBER_INTERFACE

#include "../utility/clonable.hpp"
#include "../utility/writable.hpp"

#include "number.decl.hpp"
#include "float.decl.hpp"
#include "operators.hpp"

namespace Ariadne {

/************ Number *********************************************************/

class NumberInterface;
class UnaryOperatorInterface;
template<class... YS> struct Aware;
template<class X> class NumberWrapper;



class NumberInterface
    : public std::enable_shared_from_this<NumberInterface>
    , public WritableInterface
    , public ClonableInterface
{
    template<class X> friend class NumberWrapper;
    friend class Handle<NumberInterface>;
  public:
    virtual ~NumberInterface() = default;
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

    virtual FloatDPBall _get(MetricTag, DoublePrecision, DoublePrecision) const = 0;
    virtual FloatDPBounds _get(BoundedTag, DoublePrecision) const = 0;
    virtual FloatDPUpperBound _get(UpperTag, DoublePrecision) const = 0;
    virtual FloatDPLowerBound _get(LowerTag, DoublePrecision) const = 0;
    virtual FloatDPApproximation _get(ApproximateTag, DoublePrecision) const = 0;

    virtual FloatMPDPBall _get(MetricTag, MultiplePrecision, DoublePrecision) const = 0;

    virtual FloatMPBall _get(MetricTag, MultiplePrecision, MultiplePrecision) const = 0;
    virtual FloatMPBounds _get(BoundedTag, MultiplePrecision) const = 0;
    virtual FloatMPUpperBound _get(UpperTag, MultiplePrecision) const = 0;
    virtual FloatMPLowerBound _get(LowerTag, MultiplePrecision) const = 0;
    virtual FloatMPApproximation _get(ApproximateTag, MultiplePrecision) const = 0;

    // FIXME: Only used in templates; should not really be needed
    FloatDPBall _get(MetricTag p, DoublePrecision pr) const;
    FloatMPBall _get(MetricTag p, MultiplePrecision pr) const;

    virtual ParadigmCode _paradigm() const = 0;
    virtual String _class_name() const = 0;

    virtual LogicalValue _equals(NumberInterface const& y) const = 0;
    virtual LogicalValue _less(NumberInterface const& y) const = 0;
};


} // namespace Ariadne

#endif /* ARIADNE_NUMBER_INTERFACE */
