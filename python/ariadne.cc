/***************************************************************************
 *            python/ariadne.cc
 *
 *  22 June 2005
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
#include "numerical_type.h"
#include "interval.h"
#include "state.h"
#include "rectangle.h"
#include "list_set.h"

#include <boost/python.hpp>

using namespace Ariadne;

// Since gmpxx operators return intermediate types which need conversion, need to explicitly provide arithmetic functions
inline Integer neg_z(const Integer& z1) { return -z1; }
inline Integer add_z(const Integer& z1, const Integer& z2) { return z1+z2; }
inline Integer sub_z(const Integer& z1, const Integer& z2) { return z1-z2; }
inline Integer mul_z(const Integer& z1, const Integer& z2) { return z1*z2; }
inline Rational div_z(const Integer& z1, const Integer& z2) { return Rational(z1)/z2; }
inline bool eq_z(const Integer& z1, const Integer& z2) { return z1==z2; }
inline bool ne_z(const Integer& z1, const Integer& z2) { return z1!=z2; }
inline bool lt_z(const Integer& z1, const Integer& z2) { return z1<z2; }
inline bool le_z(const Integer& z1, const Integer& z2) { return z1<=z2; }

BOOST_PYTHON_MODULE(ariadne)
{
  using boost::python::class_;
  using boost::python::init;
  using boost::python::self;
  using boost::python::return_value_policy;
  using boost::python::copy_const_reference;
  using boost::python::def;

}
