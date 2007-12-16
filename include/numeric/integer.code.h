/***************************************************************************
 *            integer.code.h
 *
 *  Copyright  2004-7  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
#include <sstream>
#include "numeric/integer.h"

namespace Ariadne {

Numeric::Integer::Integer(const std::string& str) {
  mpz_init(_value); std::stringstream ss(str); ss>>*this; }

std::ostream& 
Numeric::operator<<(std::ostream& os, const Integer& n)
{
  return os << n._value;
}

std::istream& 
Numeric::operator>>(std::istream& is, Integer& n) 
{
  return is >> n._value;
}


} // namespace Ariadne




