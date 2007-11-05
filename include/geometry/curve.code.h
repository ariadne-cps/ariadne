/***************************************************************************
 *            curve.code.h
 *
 *  Copyright  2007  Pieter Collins
 *  pieter.collins@cwi.nl
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
 
#include <cassert>

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "geometry/point.h"


namespace Ariadne {


template<class R>
Geometry::Curve<R>::~Curve() 
{
  delete this->_function_ptr;
}

template<class R>
Geometry::Curve<R>::Curve(const Function::DifferentiableFunctionInterface<R>& f) 
  : _function_ptr(f.clone())
{
  assert(this->_function_ptr->argument_size()==1);
}

template<class R>
Geometry::Curve<R>::Curve(const Curve<R>& c) 
  : _function_ptr(c._function_ptr->clone())
{
}

template<class R>
Geometry::Curve<R>* 
Geometry::Curve<R>::clone() const 
{
  return new Curve<R>(*this);
}

template<class R>
dimension_type 
Geometry::Curve<R>::dimension() const 
{
  return this->_function_ptr->result_size();
}

template<class R>
smoothness_type 
Geometry::Curve<R>::smoothness() const 
{
  return this->_function_ptr->smoothness();
}


template<class R>
Geometry::Point< typename Geometry::Curve<R>::A > 
Geometry::Curve<R>::value(const A& s) const 
{
  LinearAlgebra::Vector<A> v(1,&s);
  return Point<A>(this->_function_ptr->evaluate(v));
}

template<class R>
LinearAlgebra::Vector< typename Geometry::Curve<R>::A > 
Geometry::Curve<R>::tangent(const A& s) const 
{
  LinearAlgebra::Vector<A> v(1,&s);
  return this->_function_ptr->jacobian(v).column(0);
}


template<class R>
std::ostream& 
Geometry::Curve<R>::write(std::ostream& os) const 
{
  return os << "Curve( function=" << *this->_function_ptr << " )";
}

 
}
