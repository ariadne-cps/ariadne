/***************************************************************************
 *            evolution_parameters.h
 *
 *  Copyright  2007-8  Davide Bresolim, Alberto Casagrande, Pieter Collins
 *
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
 
/*! \file evolution_parameters.h
 *  \brief Parameters for controlling the accuracy of evaluation methods.
 */

#ifndef ARIADNE_EVOLUTION_PARAMETERS_H
#define ARIADNE_EVOLUTION_PARAMETERS_H

#include <cstddef>
#include <boost/smart_ptr.hpp>


namespace Ariadne {
  

//! \brief Parameters for controlling the accuracy of evolution methods.
//!
//! \internal This class is an abomination. Parameter specification should be
//! delegated to the class implementing the method that needs it. Currently
//! used to present a semi-stable interface to users.
struct EvolutionParameters {
  typedef int IntType;
  typedef double TimeType;
  typedef double RealType;
  
  //! \brief Default constructer gives reasonable values. 
  EvolutionParameters();

  //! \brief The time after which an evolver may approximate computed sets on a grid,
  //! in order to use previously cached results for the grid. Increasing this 
  //! parameter may improve the accuracy of the computations.  
  //! If there is recurrence in the system, then this parameter should be set to 
  //! the average recurrence time, if known. Used for discrete-time computation.
  uint lock_to_grid_steps;
    
  //! \brief A suggested minimum step size for integration. 
  //! This value may be ignored if an integration step cannot be performed without reducing the step size below this value. 
  TimeType minimum_step_size;

  //! \brief The maximum allowable step size for integration. 
  //! Decreasing this value increases the accuracy of the computation. 
  TimeType maximum_step_size;

  //! \brief A suggested minimum radius of a basic set after a subdivision (not a strict bound). 
  RealType minimum_enclosure_radius;

  //! \brief The maximum allowable radius of a basic set during integration. 
  //! Decreasing this value increases the accuracy of the computation of an over-approximation. 
  RealType maximum_enclosure_radius;

  //! \brief The time after which an VectorFieldEvolver may approximate computed sets on a grid,
  //! in order to use previously cached integration results for the grid. Increasing this 
  //! parameter improves the accuracy of the computations. Setting this parameter too
  //! low usually results in meaningless computations. A typical system trajectory 
  //! should move at least four times the grid size between locking to the grid. <br>
  //! For forced oscillators, this parameter should be set to the forcing time, 
  //! or a multiple or fraction thereof. Used for continuous-time computation.
  //!
  TimeType lock_to_grid_time;

  //! \brief Set the length of the approximation grid. 
  //! Decreasing this value increases the accuracy of the computation. 
  RealType grid_length;

  //! \brief Set the depth used for approximation on a grid
  //! Increasing this value increases the accuracy of the computation. 
  IntType grid_depth;

  //! \brief Set the size of the region used for computation. 
  //! Increasing this value reduces the risk of error due to missing orbits which leave the bounding domain. 
  RealType bounding_domain_size;

};

inline
EvolutionParameters::EvolutionParameters() 
  : lock_to_grid_steps(1),
    minimum_step_size(0.0),
    maximum_step_size(1.0),
    minimum_enclosure_radius(0.0),
    maximum_enclosure_radius(0.5),
    lock_to_grid_time(1.0),
    grid_length(0.125),
    grid_depth(6),
    bounding_domain_size(8.0)
{
}


inline
std::ostream& 
operator<<(std::ostream& os, const EvolutionParameters& p) 
{
  os << "EvolutionParameters"
     << "(\n  lock_to_grid_steps=" << p.lock_to_grid_steps

     << ",\n  minimum_step_size=" << p.minimum_step_size
     << ",\n  maximum_step_size=" << p.maximum_step_size
     << ",\n  lock_to_grid_time=" << p.lock_to_grid_time

     << ",\n  minimum_enclosure_radius=" << p.minimum_enclosure_radius
     << ",\n  maximum_enclosure_radius=" << p.maximum_enclosure_radius

     << ",\n  grid_length=" << p.grid_length

     << ",\n  bounding_domain_size=" << p.bounding_domain_size
     << "\n)\n";
  return os;
}

} //!namespace Ariadne

#endif //!ARIADNE_EVOLUTION_PARAMETERS_H
