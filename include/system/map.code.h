/***************************************************************************
 *            map.code.h
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it,  Pieter.Collins@cwi.nl
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
 
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "function/multi_index.h"
#include "function/taylor_derivative.h"
#include "function/function_interface.h"

#include "geometry/point.h"
#include "geometry/box.h"

#include "system/map.h"

#include "output/logging.h"

namespace Ariadne {

template<class R>
System::Map<R>::Map(const Function::FunctionInterface<R>& f)
  : _function_ptr(f.clone())
{
}



template<class R>
const Function::FunctionInterface<R>&
System::Map<R>::function() const
{
  return *this->_function_ptr;
}

template<class R>
Geometry::EuclideanSpace
System::Map<R>::state_space() const
{
  ARIADNE_ASSERT(this->_function_ptr->argument_size()==this->_function_ptr->result_size());
  return Geometry::EuclideanSpace(this->_function_ptr->result_size());
}

template<class R>
dimension_type 
System::Map<R>::dimension() const
{
  ARIADNE_ASSERT(this->_function_ptr->argument_size()==this->_function_ptr->result_size());
  return this->_function_ptr->result_size();
}

template<class R>
smoothness_type 
System::Map<R>::smoothness() const
{
  return this->_function_ptr->smoothness();
}

template<class R>
dimension_type 
System::Map<R>::result_dimension() const
{
  return this->_function_ptr->result_size();
}

template<class R>
dimension_type 
System::Map<R>::argument_dimension() const
{
  return this->_function_ptr->argument_size();
}





template<class R>
Geometry::Point<typename System::Map<R>::A>
System::Map<R>::operator() (const Geometry::Point<A>& x) const 
{
  ARIADNE_LOG(8,"Map::operator() (Point x) with x="<<x);
  return this->image(x);
}

template<class R>
Geometry::Point<typename System::Map<R>::A>
System::Map<R>::image(const Geometry::Point<A>& x) const 
{
  ARIADNE_LOG(8,"Map::image(Point x) with x="<<x);
  return Geometry::Point<A>(this->_function_ptr->evaluate(x.position_vector()));
}


template<class R>
LinearAlgebra::Matrix<typename System::Map<R>::A>
System::Map<R>::jacobian(const Geometry::Point<A>& x) const 
{
  return this->_function_ptr->jacobian(x.position_vector());
}

template<class R>
Function::TaylorDerivative<typename System::Map<R>::A>
System::Map<R>::derivative(const Geometry::Point<A>& x, const smoothness_type& s) const 
{
  return this->_function_ptr->derivative(x.position_vector(),s);
}

template<class R>
std::ostream&
System::Map<R>::write(std::ostream& os) const 
{
  return os << "Map( \nfunction="<<*this->_function_ptr<<"\n)";
}


}
