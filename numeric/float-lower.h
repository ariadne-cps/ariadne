/***************************************************************************
 *            float-lower.h
 *
 *  Copyright 2008-14  Pieter Collins
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

/*! \file float-lower.h
 *  \brief Floating-point lower bounds for a number.
 */

#ifndef ARIADNE_FLOAT_LOWER_H
#define ARIADNE_FLOAT_LOWER_H

#include <iostream>
#include <cassert>

#include "utility/declarations.h"

#include "utility/tribool.h"
#include "numeric/float-raw.h"

namespace Ariadne {

//! \ingroup NumericModule
//! \brief Floating-point lower bounds for real numbers.
class LowerFloat {
  public:
    typedef Lower Paradigm;
    typedef LowerFloat NumericType;
  public:
    //! \brief Default constructor yields a lower bound of \a 0.
    LowerFloat() : l(0.0) { }
    //! \brief Convert from a builtin integer.
    template<class N, EnableIf<IsIntegral<N>> = dummy> LowerFloat(N n) : l(n) { }
    //! \brief Convert from a builtin floating-point value.
    template<class X, EnableIf<IsFloatingPoint<X>> = dummy> explicit LowerFloat(X x) : l(x) { }
    //! \brief Explicitly construct from a raw floating-point value.
    explicit LowerFloat(Float x) : l(x) { }
    //! \brief Convert from floating-point bounds on a number.
    inline LowerFloat(const ValidatedFloat& x);
    //! \brief Convert from a floating-point number with an exact representation.
    LowerFloat(const ExactFloat& x);
    //! \brief Construct from a generic number.
    explicit LowerFloat(const Number<Lower>& x);
    //! \brief Assign from a generic number FIXME: Find good overloads
    // LowerFloat& operator=(Number<Lower> const& x);
    //! \brief Convert to generic number type.
    operator Number<Lower> () const;
    //! \brief Convert from a real number.
    explicit LowerFloat(const Real& x);
    explicit LowerFloat(const Rational& x);
    explicit LowerFloat(const Integer& x);
    //! \brief Explicitly convert to the raw floating-point value.
    explicit operator Float const& () const { return l; }
    //! \brief Get the raw value.
    Float const& raw() const { return l; }
    Float& raw() { return l; }
    //! \brief Get the value to double-precision.
    double get_d() const { return l.get_d(); }
    friend LowerFloat operator+(LowerFloat);
    friend UpperFloat operator-(LowerFloat);
    friend LowerFloat operator+(LowerFloat, LowerFloat);
    friend LowerFloat operator-(LowerFloat, UpperFloat);
    friend UpperFloat operator-(UpperFloat, LowerFloat);
    friend LowerFloat operator*(LowerFloat, LowerFloat);
    friend LowerFloat operator/(LowerFloat, UpperFloat);
    friend UpperFloat operator/(UpperFloat, LowerFloat);
    friend UpperFloat rec(LowerFloat);
    friend LowerFloat max(LowerFloat, LowerFloat);
    friend LowerFloat min(LowerFloat, LowerFloat);
    friend OutputStream& operator<<(OutputStream& os, LowerFloat);

    friend ApproximateFloat operator+(ApproximateFloat, ApproximateFloat);
    friend ApproximateFloat operator-(ApproximateFloat, ApproximateFloat);
    friend ApproximateFloat operator*(ApproximateFloat, ApproximateFloat);
    friend ApproximateFloat operator/(ApproximateFloat, ApproximateFloat);
  private:
    Float l;
};

} // namespace Ariadne

#endif
