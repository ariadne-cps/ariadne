/***************************************************************************
 *            float-validated.cc
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

#include "utility/standard.h"

#include <iostream>
#include <iomanip>
#include <cassert>
#include "utility/container.h"



#include "config.h"
#include "utility/typedefs.h"
#include "utility/macros.h"
#include "utility/exceptions.h"
#include "numeric/integer.h"
#include "numeric/decimal.h"
#include "numeric/dyadic.h"
#include "numeric/rational.h"
#include "numeric/real.h"
#include "numeric/number.h"
#include "numeric/float.h"
#include "numeric/float-exact.h"
#include "numeric/float-validated.h"
#include "numeric/float-approximate.h"


namespace Ariadne {

UpperFloat::UpperFloat(ExactFloat const& x) :  UpperFloat(x.raw()) {
}

UpperFloat::UpperFloat(Number<Upper> const& x) {
    ARIADNE_NOT_IMPLEMENTED;
}

UpperFloat operator+(UpperFloat x)
{
    return UpperFloat(+x.raw());
}

UpperFloat operator-(LowerFloat x)
{
    return UpperFloat(-x.raw());
}

UpperFloat operator+(UpperFloat x1, UpperFloat x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    Float x1u=x1.raw();
    Float x2u=x2.raw();
    set_rounding_mode(upward);
    Float ru=x1u+x2u;
    set_rounding_mode(rnd);
    return UpperFloat(ru);
}

UpperFloat operator-(UpperFloat x1, LowerFloat x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    Float x1u=x1.raw();
    Float x2l=x2.raw();
    set_rounding_mode(upward);
    Float ru=x1u-x2l;
    set_rounding_mode(rnd);
    return UpperFloat(ru);
}

OutputStream& operator<<(OutputStream& os, UpperFloat x) {
    rounding_mode_t rnd=get_rounding_mode();
    set_rounding_upward();
    os << std::showpoint << std::setprecision(ValidatedFloat::output_precision) << x.raw();
    set_rounding_mode(rnd);
    return os;
}

UpperFloat operator*(UpperFloat x1, UpperFloat x2) {
    assert(x1.raw()>=0.0 && x2.raw() >= 0.0);
    return UpperFloat(x1.raw()*x2.raw());
}

UpperFloat operator/(UpperFloat x1, LowerFloat x2) {
    assert(x1.raw()>=0.0 && x2.raw() > 0.0);
    return UpperFloat(x1.raw()/x2.raw());
}

UpperFloat pow(UpperFloat x, Nat n) {
    assert(x.raw()>=0.0);
    return UpperFloat(pow_up(x.raw(),n));
}

UpperFloat abs(UpperFloat x) {
    return UpperFloat(abs(Float(x)));
}

UpperFloat half(UpperFloat x) {
    return UpperFloat(half(Float(x)));
}

UpperFloat& operator+=(UpperFloat& x1, UpperFloat x2) {
     return x1=x1+x2;
}

UpperFloat& operator*=(UpperFloat& x1, UpperFloat x2) {
     return x1=x1*x2;
}

UpperFloat rec(LowerFloat x)
{
    rounding_mode_t rnd=get_rounding_mode();
    Float xl=x.raw();
    set_rounding_mode(upward);
    Float ru=1/xl;
    set_rounding_mode(rnd);
    return UpperFloat(ru);
}


UpperFloat min(UpperFloat x1, UpperFloat x2) {
    return UpperFloat(min(x1.raw(),x2.raw()));
}

UpperFloat max(UpperFloat x1, UpperFloat x2) {
    return UpperFloat(max(x1.raw(),x2.raw()));
}




PositiveUpperFloat operator+(PositiveUpperFloat x1, PositiveUpperFloat x2) {
    return PositiveUpperFloat(add_up(x1.raw(),x2.raw()));
}

PositiveUpperFloat operator-(PositiveUpperFloat x1, LowerFloat x2) {
    ARIADNE_PRECONDITION(x1.raw()>=x2.raw());
    return PositiveUpperFloat(sub_up(x1.raw(),x2.raw()));
}

PositiveUpperFloat operator*(PositiveUpperFloat x1, PositiveUpperFloat x2) {
    return PositiveUpperFloat(mul_up(x1.raw(),x2.raw()));
}

PositiveUpperFloat operator/(PositiveUpperFloat x1, LowerFloat x2) {
    ARIADNE_PRECONDITION(x2.raw()>=0);
    return PositiveUpperFloat(div_up(x1.raw(),x2.raw()));
}

PositiveUpperFloat pow(PositiveUpperFloat x, Nat m) {
    return PositiveUpperFloat(pow_up(x.raw(),m));
}

PositiveUpperFloat half(PositiveUpperFloat x) {
    return PositiveUpperFloat(half(x.raw()));
}

} // namespace Ariadne

