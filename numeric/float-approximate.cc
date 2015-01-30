/***************************************************************************
 *            float-approximate.cc
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
#include <limits>



#include "config.h"

#include "utility/macros.h"
#include "numeric/dyadic.h"
#include "numeric/decimal.h"
#include "numeric/rational.h"

#include "numeric/float.h"
#include "numeric/float-exact.h"
#include "numeric/float-validated.h"
#include "numeric/float-approximate.h"


namespace Ariadne {


Nat ApproximateFloat::output_precision = 6;
Nat ExactFloat::output_precision = 18;

const ExactFloat infty = ExactFloat(inf);


//ExactFloat inf = ExactFloat(std::numeric_limits< double >::infinity());
ApproximateFloat::ApproximateFloat(Dyadic const& b) : ApproximateFloat(b.operator Rational()) { }
ApproximateFloat::ApproximateFloat(Decimal const& d) : ApproximateFloat(d.operator Rational()) { }

ApproximateFloat::ApproximateFloat(ExactFloat const& x) : ApproximateFloat(x.raw()) { }
ApproximateFloat::ApproximateFloat(ValidatedFloat const& x) : ApproximateFloat(half_exact(add_near(x.lower_raw(),x.upper_raw()))) { }
ApproximateFloat::ApproximateFloat(UpperFloat const& x) : ApproximateFloat(x.raw()) { }
ApproximateFloat::ApproximateFloat(LowerFloat const& x) : ApproximateFloat(x.raw()) { }

ApproximateFloat::ApproximateFloat(Number<Approximate> const& x) { ARIADNE_NOT_IMPLEMENTED; }

ExactFloat::operator Rational() const {
    return Rational(this->get_d());
}
ApproximateFloat::ApproximateFloat(Rational const& q) : ApproximateFloat(q.get_d()) {
}

} // namespace Ariadne
