/***************************************************************************
 *            vector_field.code.h
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
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
 
#include "vector_field.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "function/taylor_derivative.h"
#include "function/function_interface.h"

#include "geometry/point.h"
#include "geometry/box.h"



namespace Ariadne {

template<class R>
System::VectorField<R>::~VectorField() 
{
}
     
template<class R>
System::VectorField<R>::VectorField(const Function::FunctionInterface<R>& f) 
  : _function_ptr(f.clone())
{
  assert(f.argument_size()==f.result_size());
}
     

template<class R>
LinearAlgebra::Vector<typename System::VectorField<R>::F> 
System::VectorField<R>::evaluate(const Geometry::Point<F>& x) const 
{
  return this->_function_ptr->evaluate(x.position_vector());
}

template<class R>
LinearAlgebra::Matrix<typename System::VectorField<R>::F> 
System::VectorField<R>::jacobian(const Geometry::Point<F>& x) const 
{
  return this->_function_ptr->jacobian(x.position_vector());
}

template<class R>
Function::TaylorDerivative<typename System::VectorField<R>::F>
System::VectorField<R>::derivative(const Geometry::Point<F>& x, const smoothness_type& s) const
{
  return this->_function_ptr->derivative(x.position_vector(),s);
}
   
template<class R>
std::ostream&
System::VectorField<R>::write(std::ostream& os) const 
{
  return os << "VectorField( function=" << *this->_function_ptr << " )";
}

    

}
