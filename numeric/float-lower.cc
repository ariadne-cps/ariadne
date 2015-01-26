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
#include "numeric/float-upper.h"
#include "numeric/float-lower.h"
#include "numeric/float-approximate.h"


namespace Ariadne {

LowerFloat::LowerFloat(ExactFloat const& x) :  LowerFloat(x.raw()) {
}

LowerFloat::LowerFloat(Number<Lower> const& x) {
    ARIADNE_NOT_IMPLEMENTED;
}

LowerFloat operator+(LowerFloat x)
{
    volatile double xl=internal_cast<volatile double&>(x.raw());
    volatile double rl=+xl;
    return LowerFloat(rl);
}

LowerFloat operator-(UpperFloat x)
{
    volatile double xu=internal_cast<volatile double&>(x.raw());
    volatile double rl=-xu;
    return LowerFloat(rl);
}

LowerFloat operator+(LowerFloat x1, LowerFloat x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double x1l=internal_cast<volatile double&>(x1.raw());
    volatile double x2l=internal_cast<volatile double&>(x2.raw());
    set_rounding_mode(downward);
    volatile double rl=x1l+x2l;
    set_rounding_mode(rnd);
    return LowerFloat(rl);
}

LowerFloat operator-(LowerFloat x1, UpperFloat x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double x1l=internal_cast<volatile double&>(x1.raw());
    volatile double x2u=internal_cast<volatile double&>(x2.raw());
    set_rounding_mode(downward);
    volatile double rl=x1l-x2u;
    set_rounding_mode(rnd);
    return LowerFloat(rl);
}

LowerFloat operator*(LowerFloat x1, LowerFloat x2)
{
    ARIADNE_PRECONDITION(x1.raw()>=0 && x2.raw()>=0);
    rounding_mode_t rnd=get_rounding_mode();
    volatile double x1l=internal_cast<volatile double&>(x1.raw());
    volatile double x2l=internal_cast<volatile double&>(x2.raw());
    set_rounding_mode(downward);
    volatile double rl=x1l*x2l;
    set_rounding_mode(rnd);
    return LowerFloat(rl);
}

OutputStream& operator<<(OutputStream& os, LowerFloat x) {
    rounding_mode_t rnd=get_rounding_mode();
    set_rounding_downward();
    os << std::showpoint << std::setprecision(ValidatedFloat::output_precision) << x.raw();
    set_rounding_mode(rnd);
    return os;
}

LowerFloat rec(UpperFloat x)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double xu=internal_cast<volatile double&>(x.raw());
    set_rounding_mode(downward);
    volatile double rl=1/xu;
    set_rounding_mode(rnd);
    return LowerFloat(rl);
}
LowerFloat min(LowerFloat x1, LowerFloat x2) {
    return LowerFloat(min(x1.raw(),x2.raw()));
}

LowerFloat max(LowerFloat x1, LowerFloat x2) {
    return LowerFloat(max(x1.raw(),x2.raw()));
}

} // namespace Ariadne

