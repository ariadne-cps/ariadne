/***************************************************************************
 *            float.cc
 *
 *  Copyright 2008-10  Pieter Collins
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

#include <iostream>
#include <iomanip>
#include <cassert>
#include <limits>

#include "config.h"

#include "macros.h"
#include "dyadic.h"
#include "decimal.h"
#include "rational.h"

#include "float.h"
#include "float-exact.h"
#include "float-validated.h"
#include "float-approximate.h"


namespace Ariadne {


uint ApproximateFloat::output_precision = 6;
uint ExactFloat::output_precision = 18;

//ExactFloat inf = ExactFloat(std::numeric_limits< double >::infinity());
ApproximateFloat::ApproximateFloat(Dyadic const& b) : ApproximateFloat(b.operator Rational()) { }
ApproximateFloat::ApproximateFloat(Decimal const& d) : ApproximateFloat(d.operator Rational()) { }

ApproximateFloat::ApproximateFloat(ExactFloat const& x) : ApproximateFloat(x.operator Float()) { }
ApproximateFloat::ApproximateFloat(ValidatedFloat const& x) : ApproximateFloat(half_exact(add_near(x.lower(),x.upper()))) { }

#ifdef HAVE_GMPXX_H
ExactFloat::operator Rational() const {
    return Rational(this->get_d());
}
ApproximateFloat::ApproximateFloat(Rational const& q) : ApproximateFloat(q.get_d()) {
}
#endif

} // namespace Ariadne


#ifdef ARIADNE_ENABLE_SERIALIZATION

#include "serialization.h"

namespace Ariadne {

void serialize(boost::archive::text_oarchive& a, Float& flt, const unsigned int v) {
    const double x=flt.get_d();
    a << x;
};

void serialize(boost::archive::text_iarchive& a, Float& flt, const unsigned int v) {
    flt=std::numeric_limits<double>::quiet_NaN();
    double x;
    a >> x;
    flt=x;
}

} // namespace Ariadne

#endif /* ARIADNE_ENABLE_SERIALIZATION */



