/***************************************************************************
 *            polynomial_map.code.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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

#include "polynomial_map.h"

#include <string>
#include <sstream>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "base/stlio.h"
#include "base/array.h"
#include "numeric/interval.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "geometry/point.h"
#include "geometry/box.h"

#include "system/exceptions.h"

namespace Ariadne {




template<class R>
Geometry::Point<typename System::PolynomialMap<R>::F>
System::PolynomialMap<R>::image(const Geometry::Point<F>& s) const 
{
  ARIADNE_CHECK_ARGUMENT_DIMENSION(*this,s,"System::PolynomialMap<R>::apply(Point<R>)");
  return Geometry::Point<F>(this->_polynomial.evaluate(s.position_vector()));
}


template<class R>
LinearAlgebra::Matrix< typename System::PolynomialMap<R>::F >
System::PolynomialMap<R>::jacobian(const Geometry::Point<F>& s) const 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}







template<class R>
std::ostream&
System::PolynomialMap<R>::write(std::ostream& os) const
{
  const PolynomialMap<R>& pm = *this;
  return os << pm._polynomial;
}


template<class R>
std::istream&
System::PolynomialMap<R>::read(std::istream& is)
{
  PolynomialMap<R>& pm=*this;
  is >> pm._polynomial;
  return is;
}


}
