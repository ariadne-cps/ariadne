/***************************************************************************
 *            evaluation_parameters.code.h
 *
 *  Copyright  2007-8  Davide Bresolin, Alberto Casagrande, Pieter Collins
 *  davide.bresolin@univr.it, casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
#include "geometry/box.h"
#include "geometry/grid.h"
#include "geometry/hybrid_denotable_set.h"
#include "geometry/hybrid_set.h"

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
    _bounding_domain_size(1),
    _verbosity(0),
    _grid(0,1),
    _hybrid_grid(),
		_hybrid_bounding_domain()
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
uint
Evaluation::EvolutionParameters<R>::verbosity() const 
{
  return this->_verbosity;
}




template<class R>
Geometry::Box<R>
Evaluation::EvolutionParameters<R>::bounding_domain(dimension_type d) const 
{
  return Geometry::Box<R>(LinearAlgebra::Vector< Numeric::Interval<R> >(d,Numeric::Interval<R>(-1,1)*this->bounding_domain_size()));
}


template<class R>
Geometry::HybridSet<R>
Evaluation::EvolutionParameters<R>::hybrid_bounding_domain(const Geometry::HybridSpace& loc) const
{
  if(this->_hybrid_bounding_domain.locations() == loc) {
    return this->_hybrid_bounding_domain;
  }
	Geometry::HybridSet<R> result;
	for(typename Geometry::HybridSpace::const_iterator loc_iter=loc.begin();
			loc_iter!=loc.end(); ++loc_iter)
	{
		result.new_location(loc_iter->discrete_state(),Geometry::RectangularSet<R>(this->bounding_domain(loc_iter->dimension())));
	}
	return result;
}


template<class R>
Geometry::Grid<R>
Evaluation::EvolutionParameters<R>::grid(dimension_type d) const 
{
  if(this->_grid.dimension() == d) {
    return this->_grid;
  }
  return Geometry::Grid<R>(LinearAlgebra::Vector<R>(d,this->_grid_length));
}


template<class R>
Geometry::FiniteGrid<R>
Evaluation::EvolutionParameters<R>::finite_grid(dimension_type d) const 
{
  return Geometry::FiniteGrid<R>(this->grid(d),this->bounding_domain(d));
}

template<class R>
Geometry::HybridGrid<R>
Evaluation::EvolutionParameters<R>::hybrid_grid(const Geometry::HybridSpace& loc) const 
{
  if(this->_hybrid_grid.locations() == loc) {
    return this->_hybrid_grid;
  }
  return Geometry::HybridGrid<R>(loc,this->_grid_length);
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
Evaluation::EvolutionParameters<R>::set_minimum_step_size(double x) 
{
  this->_minimum_step_size=x;
}

template<class R>
void
Evaluation::EvolutionParameters<R>::set_maximum_step_size(Numeric::Rational x)  
{
  this->_maximum_step_size=x;
}

template<class R>
void
Evaluation::EvolutionParameters<R>::set_maximum_step_size(double x)  
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
Evaluation::EvolutionParameters<R>::set_lock_to_grid_time(double x)  
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
Evaluation::EvolutionParameters<R>::set_minimum_basic_set_radius(double x) 
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
Evaluation::EvolutionParameters<R>::set_maximum_basic_set_radius(double x)  
{
  this->_maximum_basic_set_radius=x;
}


template<class R>
void
Evaluation::EvolutionParameters<R>::set_grid_length(double x)  
{
  this->_grid_length=x;
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
Evaluation::EvolutionParameters<R>::set_argument_grid_length(double x)  
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
Evaluation::EvolutionParameters<R>::set_result_grid_length(double x) 
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
void
Evaluation::EvolutionParameters<R>::set_bounding_domain_size(double x)  
{
  this->_bounding_domain_size=x;
}



template<class R>
void
Evaluation::EvolutionParameters<R>::set_verbosity(uint v)  
{
  this->_verbosity=v;
}



template<class R>
void
  Evaluation::EvolutionParameters<R>::set_hybrid_grid(Geometry::HybridGrid<R> hgrid)  
{
  this->_hybrid_grid=hgrid;
}


template<class R>
void
  Evaluation::EvolutionParameters<R>::set_hybrid_bounding_domain(Geometry::HybridSet<R> hdom)  
{
  this->_hybrid_bounding_domain=hdom;
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
		 << ",\n  grid=" << this->_grid
		 << ",\n  hybrid_grid=" << this->_hybrid_grid

     << ",\n  bounding_domain_size=" << this->_bounding_domain_size
		 << ",\n  hybrid_bounding_domain=" << this->_hybrid_bounding_domain
     << ",\n  verbosity=" << this->_verbosity
     << "\n)\n";
  return os;
}


}
