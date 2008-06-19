/***************************************************************************
 *            solver.inline.h
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
 

#include <exception>
#include <stdexcept>
#include <string>

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "geometry/point.h"
#include "geometry/box.h"

#include "function/taylor_derivative.h"
#include "system/map.h"

namespace Ariadne {


template<class R> 
inline
DifferenceFunction<R>::DifferenceFunction(const FunctionInterface<R>& f)
  : _base(f)
{ 
  if(f.argument_size()!=f.result_size()) { 
    throw IncompatibleDimensions("DifferenceDifferenceFunction(Map f): The argument dimension must equal the result dimension"); 
  } 
}


template<class R> 
inline
DifferenceFunction<R>* 
DifferenceFunction<R>::clone() const 
{ 
  return new DifferenceFunction<R>(this->_base); 
}


template<class R> 
inline
smoothness_type 
DifferenceFunction<R>::smoothness() const 
{ 
  return _base.smoothness(); 
}


template<class R> 
inline
size_type 
DifferenceFunction<R>::argument_size() const 
{ 
  return _base.argument_size(); 
}

template<class R> 
inline
size_type 
DifferenceFunction<R>::result_size() const 
{ 
  return _base.result_size(); 
}


template<class R> 
inline
Vector<typename DifferenceFunction<R>::F> 
DifferenceFunction<R>::evaluate(const Vector<F>& p) const 
{
  return _base.evaluate(p)-p; 
}


template<class R> 
inline
Matrix<typename DifferenceFunction<R>::F>
DifferenceFunction<R>::jacobian(const Vector<F>& p) const 
{
  Matrix<F> d=_base.jacobian(p);
  Matrix<F> i=Matrix< Interval<R> >::identity(this->result_size());
  return d-i; 
}

template<class R> 
inline
TaylorDerivative<typename DifferenceFunction<R>::F>
DifferenceFunction<R>::derivative(const Vector<F>& p, const smoothness_type& s) const 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R> 
inline
std::string 
DifferenceFunction<R>::name() const 
{ 
  return "DifferenceFunction"; 
}

template<class R> 
inline
std::ostream&
DifferenceFunction<R>::write(std::ostream& os) const 
{ 
  return os << "DifferenceFunction"; 
}






template<class R> 
inline
SolverInterface<R>::~SolverInterface() 
{ 
}


template<class R>
inline
SolverInterface<R>::SolverInterface(R max_error, uint max_steps)
  : _max_error(max_error), _max_steps(max_steps) 
{
}


template<class R> 
inline
const R& 
SolverInterface<R>::maximum_error() const 
{
  return this->_max_error; 
}


template<class R> 
inline
const uint& 
SolverInterface<R>::maximum_number_of_steps() const 
{ 
  return this->_max_steps; 
}


template<class R> 
inline
void 
SolverInterface<R>::set_maximum_error(R max_error)
{ 
  this->_max_error=max_error; 
}


template<class R> 
inline
void 
SolverInterface<R>::set_maximum_number_of_steps(uint max_steps) 
{ 
  this->_max_steps=max_steps; 
}


template<class R> 
inline
Point<typename SolverInterface<R>::I> 
SolverInterface<R>::fixed_point(const Map<R>& f,const Point<I>& pt) 
{
  return this->solve(DifferenceFunction<R>(f.function()),pt); 
}



} // namespace Ariadne
