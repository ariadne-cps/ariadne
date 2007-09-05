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

#include "../base/stlio.h"
#include "../base/array.h"
#include "../numeric/interval.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../geometry/point.h"
#include "../geometry/rectangle.h"

#include "../system/exceptions.h"

namespace Ariadne {

template<class R>
void
System::PolynomialMap<R>::_set_argument_dimension(const dimension_type& n) 
{
  this->_argument_dimension=n;
  for(size_type i=0; i!=this->_components.size(); ++i) {
    this->_components[i]=Function::Polynomial<R>(n);
  }
}

template<class R>
dimension_type
System::PolynomialMap<R>::_compute_maximum_component_dimension() const
{
  size_type result=0;
  for(size_type i=0; i!=this->_components.size(); ++i) {
    result=std::max(this->_components[i].argument_size(),result);
  }
  return result;
}






template<class R>
Geometry::Point<typename System::PolynomialMap<R>::F>
System::PolynomialMap<R>::image(const Geometry::Point<F>& s) const 
{
  ARIADNE_CHECK_ARGUMENT_DIMENSION(*this,s,"System::PolynomialMap<R>::apply(Point<R>)");
  LinearAlgebra::Vector<F> v=s.position_vector();
  Geometry::Point<F> result(this->result_dimension());
  for(size_type i=0; i!=this->result_dimension(); ++i) {
    result[i] = _components[i].evaluate(v);
  }
  return result;
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
  return os << pm._components;
}


template<class R>
std::istream&
System::PolynomialMap<R>::read(std::istream& is)
{
  PolynomialMap<R>& pm=*this;
  std::vector< Function::Polynomial<R> > vec;
  is >> vec;
  pm._components=array< Function::Polynomial<R> >(vec.begin(),vec.end());
  pm._set_argument_dimension(pm._compute_maximum_component_dimension());
  return is;
}


}
