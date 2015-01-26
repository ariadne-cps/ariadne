/***************************************************************************
 *            float-upper.h
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

/*! \file float-upper.h
 *  \brief Floating-point upper bounds for a number.
 */

#ifndef ARIADNE_FLOAT_UPPER_H
#define ARIADNE_FLOAT_UPPER_H

#include <iostream>
#include <cassert>

#include "utility/declarations.h"

#include "utility/tribool.h"
#include "numeric/rounding.h"
#include "numeric/float.decl.h"
#include "numeric/float-raw.h"

namespace Ariadne {

//! \ingroup NumericModule
//! \brief Floating-point upper bounds for real numbers.
class UpperFloat {
  public:
    typedef Upper Paradigm;
    typedef UpperFloat NumericType;
  public:
    //! \brief Default constructor yields an upper bound of \a 0.
    UpperFloat() : u(0.0) { }
    //! \brief Convert from a builtin integer.
    template<class N, EnableIf<IsIntegral<N>> = dummy> UpperFloat(N n) : u(n) { }
    //! \brief Convert from a builtin floating-point value.
    template<class X, EnableIf<IsFloatingPoint<X>> = dummy> explicit UpperFloat(X x) : u(x) { }
    //! \brief Explicitly construct from a raw floating-point value.
    explicit UpperFloat(Float x) : u(x) { }
    //! \brief Convert from floating-point bounds on a number.
    inline UpperFloat(const ValidatedFloat& x);
    //! \brief Convert from a real number.
    explicit UpperFloat(const Real& x);
    explicit UpperFloat(const Rational& x);
    explicit UpperFloat(const Integer& x);
    //! \brief Construct from a generic number.
    explicit UpperFloat(const Number<Upper>& x);
    //! \brief Assign from a generic number FIXME: Find good overloads
    // UpperFloat& operator=(Number<Upper> const& x);
    //! \brief Convert to generic number type.
    operator Number<Upper> () const;
    //! \brief Convert from a floating-point number with an exact representation.
    UpperFloat(const ExactFloat& x);
    //! \brief Explicitly convert to the raw floating-point value.
    explicit operator Float const& () const { return u; }
    //! \brief Get the raw value.
    Float const& raw() const { return u; }
    Float& raw() { return u; }
    //! \brief Get the value to double-precision.
    double get_d() const { return u.get_d(); }
    friend UpperFloat operator+(UpperFloat);
    friend LowerFloat operator-(UpperFloat);
    friend UpperFloat operator+(UpperFloat, UpperFloat);
    friend UpperFloat operator-(UpperFloat, LowerFloat);
    friend LowerFloat operator-(LowerFloat, UpperFloat);
    friend UpperFloat operator*(UpperFloat, UpperFloat);
    friend UpperFloat operator/(UpperFloat, LowerFloat);
    friend LowerFloat operator/(LowerFloat, UpperFloat);
    friend UpperFloat& operator+=(UpperFloat&, UpperFloat);
    friend UpperFloat& operator*=(UpperFloat&, UpperFloat);
    friend UpperFloat& operator/=(UpperFloat&, Nat);
    friend LowerFloat rec(UpperFloat);
    friend UpperFloat pow(UpperFloat, Nat);
    friend UpperFloat half(UpperFloat);
    friend UpperFloat abs(UpperFloat);
    friend UpperFloat max(UpperFloat, UpperFloat);
    friend UpperFloat min(UpperFloat, UpperFloat);
    friend OutputStream& operator<<(OutputStream& os, UpperFloat);

    friend ApproximateFloat operator+(ApproximateFloat, ApproximateFloat);
    friend ApproximateFloat operator-(ApproximateFloat, ApproximateFloat);
    friend ApproximateFloat operator*(ApproximateFloat, ApproximateFloat);
    friend ApproximateFloat operator/(ApproximateFloat, ApproximateFloat);
  private:
    Float u;
};


} // namespace Ariadne

#endif
