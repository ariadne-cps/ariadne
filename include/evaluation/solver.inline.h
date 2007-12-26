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
#include "system/vector_field.h"
#include "system/map.h"

namespace Ariadne {
namespace System {

template<class R> 
inline
DifferenceMap<R>::DifferenceMap(const MapInterface<R>& f)
  : _base(f)
{ 
  if(f.argument_dimension()!=f.result_dimension()) { 
    throw Geometry::IncompatibleDimensions("DifferenceMap::DifferenceMap(Map f): The argument dimension must equal the result dimension"); 
  } 
}


template<class R> 
inline
DifferenceMap<R>* 
DifferenceMap<R>::clone() const 
{ 
  return new DifferenceMap<R>(this->_base); 
}


template<class R> 
inline
smoothness_type 
DifferenceMap<R>::smoothness() const 
{ 
  return _base.smoothness(); 
}


template<class R> 
inline
dimension_type 
DifferenceMap<R>::dimension() const 
{ 
  return _base.argument_dimension(); 
}


template<class R> 
inline
LinearAlgebra::Vector<typename DifferenceMap<R>::F> 
DifferenceMap<R>::image(const Geometry::Point<F>& p) const 
{
  return _base.image(p)-p; 
}


template<class R> 
inline
LinearAlgebra::Matrix< Numeric::Interval<R> > 
DifferenceMap<R>::jacobian(const Geometry::Point<F>& p) const 
{
  LinearAlgebra::Matrix<F> d=_base.jacobian(p);
  LinearAlgebra::Matrix<F> i=LinearAlgebra::Matrix< Numeric::Interval<R> >::identity(this->dimension());
  return d-i; 
}


template<class R> 
inline
std::string 
DifferenceMap<R>::name() const 
{ 
  return "DifferenceMap"; 
}



} namespace Evaluation {



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
Geometry::Point<typename SolverInterface<R>::I> 
SolverInterface<R>::fixed_point(const System::MapInterface<R>& f,const Geometry::Point<I>& pt) 
{
  return this->solve(System::DifferenceMap<R>(f),pt); 
}



}}
