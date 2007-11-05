/***************************************************************************
 *            evaluation_parameters.code.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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
 
#include "../linear_algebra/vector.h"
#include "../geometry/rectangle.h"
#include "../geometry/grid.h"

namespace Ariadne {

template<class R>
Evaluation::EvolutionParameters<R>::EvolutionParameters() 
  : _maximum_number_of_steps(255),
    _lock_to_grid_steps(1),
    _minimum_step_size(0),
    _maximum_step_size(1),
    _lock_to_grid_time(1),
    _minimum_basic_set_radius(0),
    _maximum_basic_set_radius(1),
    _grid_length(1),
    _argument_grid_length(1),
    _result_grid_length(1),
    _bounding_domain_size(1)
{
}


template<class R>
Evaluation::EvolutionParameters<R>*
Evaluation::EvolutionParameters<R>::clone() const
{
  return new EvolutionParameters<R>(*this);
}


template<class R>
size_type
Evaluation::EvolutionParameters<R>::maximum_number_of_steps() const
{
  return this->_maximum_number_of_steps;
}


template<class R>
size_type
Evaluation::EvolutionParameters<R>::lock_to_grid_steps() const
{
  return this->_lock_to_grid_steps;
}


template<class R>
time_type
Evaluation::EvolutionParameters<R>::minimum_step_size() const
{
  return this->_minimum_step_size;
}


template<class R>
time_type
Evaluation::EvolutionParameters<R>::maximum_step_size() const 
{
  return this->_maximum_step_size;
}


template<class R>
time_type
Evaluation::EvolutionParameters<R>::lock_to_grid_time() const 
{
  return this->_lock_to_grid_time;
}


template<class R>
R
Evaluation::EvolutionParameters<R>::minimum_basic_set_radius() const
{
  return this->_minimum_basic_set_radius;
}


template<class R>
R
Evaluation::EvolutionParameters<R>::maximum_basic_set_radius() const 
{
  return this->_maximum_basic_set_radius;
}


template<class R>
R
Evaluation::EvolutionParameters<R>::grid_length() const 
{
  return this->_grid_length;
}


template<class R>
R
Evaluation::EvolutionParameters<R>::argument_grid_length() const 
{
  return this->_argument_grid_length;
}

template<class R>
R
Evaluation::EvolutionParameters<R>::result_grid_length() const 
{
  return this->_result_grid_length;
}

template<class R>
R
Evaluation::EvolutionParameters<R>::bounding_domain_size() const 
{
  return this->_bounding_domain_size;
}




template<class R>
Geometry::Rectangle<R>
Evaluation::EvolutionParameters<R>::bounding_box(dimension_type d) const 
{
  return Geometry::Rectangle<R>(LinearAlgebra::Vector< Numeric::Interval<R> >(d,Numeric::Interval<R>(-1,1)*this->bounding_domain_size()));
}


template<class R>
Geometry::Grid<R>
Evaluation::EvolutionParameters<R>::grid(dimension_type d) const 
{
  return Geometry::Grid<R>(LinearAlgebra::Vector<R>(d,this->_grid_length));
}


template<class R>
Geometry::FiniteGrid<R>
Evaluation::EvolutionParameters<R>::finite_grid(dimension_type d) const 
{
  return Geometry::FiniteGrid<R>(this->grid(d),this->bounding_box(d));
}




template<class R>
void
Evaluation::EvolutionParameters<R>::set_maximum_number_of_steps(size_type x) 
{
  this->_maximum_number_of_steps=x;
}


template<class R>
void
Evaluation::EvolutionParameters<R>::set_lock_to_grid_steps(size_type x) 
{
  this->_lock_to_grid_steps=x;
}


template<class R>
void
Evaluation::EvolutionParameters<R>::set_minimum_step_size(time_type x) 
{
  this->_minimum_step_size=x;
}


template<class R>
void
Evaluation::EvolutionParameters<R>::set_maximum_step_size(time_type x)  
{
  this->_maximum_step_size=x;
}


template<class R>
void
Evaluation::EvolutionParameters<R>::set_lock_to_grid_time(time_type x)  
{
  this->_lock_to_grid_time=x;
}


template<class R>
void
Evaluation::EvolutionParameters<R>::set_minimum_basic_set_radius(R x) 
{
  this->_minimum_basic_set_radius=x;
}


template<class R>
void
Evaluation::EvolutionParameters<R>::set_maximum_basic_set_radius(R x)  
{
  this->_maximum_basic_set_radius=x;
}


template<class R>
void
Evaluation::EvolutionParameters<R>::set_grid_length(R x)  
{
  this->_grid_length=x;
}


template<class R>
void
Evaluation::EvolutionParameters<R>::set_argument_grid_length(R x)  
{
  this->_argument_grid_length=x;
}

template<class R>
void
Evaluation::EvolutionParameters<R>::set_result_grid_length(R x) 
{
  this->_result_grid_length=x;
}

template<class R>
void
Evaluation::EvolutionParameters<R>::set_bounding_domain_size(R x)  
{
  this->_bounding_domain_size=x;
}


template<class R>
std::ostream&
Evaluation::EvolutionParameters<R>::write(std::ostream& os) const
{
  os << "EvolutionParameters"
     << "(\n  maximum_number_of_steps=" << this->_maximum_number_of_steps
     << ",\n  lock_to_grid_steps=" << this->_lock_to_grid_steps

     << ",\n  minimum_step_size=" << this->_minimum_step_size
     << ",\n  maximum_step_size=" << this->_maximum_step_size
     << ",\n  lock_to_grid_time=" << this->_lock_to_grid_time

     << ",\n  minimum_basic_set_radius=" << this->_minimum_basic_set_radius
     << ",\n  maximum_basic_set_radius=" << this->_maximum_basic_set_radius

     << ",\n  grid_length=" << this->_grid_length
     << ",\n  argument_grid_length=" << this->_argument_grid_length
     << ",\n  result_grid_length=" << this->_result_grid_length

     << ",\n  bounding_domain_size=" << this->_bounding_domain_size
     << "\n)\n";
  return os;
}


}
