/***************************************************************************
 *            function_map.code.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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
 
#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../function/function_interface.h"
#include "../function/interpreted_function.h"
#include "../function/multi_index.h"

#include "../geometry/point.h"
#include "../geometry/rectangle.h"

#include "../system/map.h"

#include "../output/logging.h"

namespace Ariadne {

template<class R>
System::FunctionMap<R>::FunctionMap(const Function::DifferentiableFunctionInterface<R>& f, const Geometry::Point<A>& param)
  : _function_ptr(f.clone()), _parameters(param)
{
  if(f.argument_size()<=param.dimension()) {
    ARIADNE_THROW(std::runtime_error,"FunctionMap(DifferentiableFunctionInterface f, Point param)",
                  "with f.argument_size()="<<f.argument_size()<<", param="<<param<<": too many parameters");
  }
}

template<class R>
System::FunctionMap<R>::FunctionMap(const Function::InterpretedFunction<R>& f, const Geometry::Point<A>& param)
  : _function_ptr(f.clone()), _parameters(param)
{
  if(f.argument_size()<=param.dimension()) {
    ARIADNE_THROW(std::runtime_error,"FunctionMap(Function f, Point param)",
                  "with f.argument_size()="<<f.argument_size()<<", param="<<param<<": too many parameters");
  }
}

template<class R>
std::string
System::FunctionMap<R>::name() const
{
  return this->_function_ptr->name();
}

template<class R>
System::FunctionMap<R>* 
System::FunctionMap<R>::clone() const
{
  return new FunctionMap<R>(*this);
}

template<class R>
smoothness_type 
System::FunctionMap<R>::smoothness() const
{
  return this->_function_ptr->smoothness();
}

template<class R>
dimension_type 
System::FunctionMap<R>::result_dimension() const
{
  return this->_function_ptr->result_size();
}

template<class R>
dimension_type 
System::FunctionMap<R>::argument_dimension() const
{
  return this->_function_ptr->argument_size()-this->_parameters.dimension();
}

template<class R>
dimension_type 
System::FunctionMap<R>::number_of_parameters() const
{
  return this->_parameters.dimension();
}


template<class R>
void 
System::FunctionMap<R>::set_parameters(const Geometry::Point<A>& param) 
{
  ARIADNE_CHECK_DIMENSION(param,this->number_of_parameters(),"FunctionMap::set_parameters(Point param)");
  this->_parameters=param;
}

template<class R>
const Geometry::Point<typename System::FunctionMap<R>::A>& 
System::FunctionMap<R>::parameters() const 
{
  return this->_parameters;
}




template<class R>
Geometry::Point<typename System::FunctionMap<R>::A>
System::FunctionMap<R>::image(const Geometry::Point<A>& x) const 
{
  ARIADNE_LOG(8,"FunctionMap::image(Point x) with x="<<x);
  LinearAlgebra::Vector<A> v=LinearAlgebra::direct_sum(x.position_vector(),this->_parameters.position_vector());
  return Geometry::Point<A>(this->_function_ptr->evaluate(v));
}


template<class R>
LinearAlgebra::Matrix<typename System::FunctionMap<R>::A>
System::FunctionMap<R>::jacobian(const Geometry::Point<A>& x) const 
{
  size_type m=this->result_dimension();
  size_type n=this->argument_dimension();
  size_type np=this->number_of_parameters();
  LinearAlgebra::Vector<A> v=direct_sum(x.position_vector(),this->_parameters.position_vector());
  LinearAlgebra::Matrix<A> J=this->_function_ptr->jacobian(v);
  // assume row major
  return LinearAlgebra::Matrix<A>(m,n,J.begin(),np,1u);
}

template<class R>
typename System::FunctionMap<R>::A
System::FunctionMap<R>::derivative(const Geometry::Point<A>& x, const size_type& i, const Function::MultiIndex& j) const 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class R>
std::ostream&
System::FunctionMap<R>::write(std::ostream& os) const 
{
  return os << "FunctionMap( \nfunction="<<*this->_function_ptr<<",\nparameters="<<this->_parameters<<"\n)";
}


}

