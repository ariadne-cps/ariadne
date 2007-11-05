/***************************************************************************
 *            evolution_parameters.h
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
 
/*! \file evolution_parameters.h
 *  \brief Parameters for controlling the accuracy of evaluation methods.
 */

#ifndef ARIADNE_EVOLUTION_PARAMETERS_H
#define ARIADNE_EVOLUTION_PARAMETERS_H

#include <boost/smart_ptr.hpp>

#include "../base/types.h"
#include "../geometry/declarations.h"

namespace Ariadne {
  namespace Evaluation {

    /*! \brief Parameters for controlling the accuracy of evolution methods. */
    template<class R>
    class EvolutionParameters {
     private:
      size_type _maximum_number_of_steps;
      size_type _lock_to_grid_steps;

      time_type _minimum_step_size;
      time_type _maximum_step_size;
      time_type _lock_to_grid_time;

      R _minimum_basic_set_radius;
      R _maximum_basic_set_radius;

      R _grid_length;
      R _argument_grid_length;
      R _result_grid_length;

      R _bounding_domain_size;
     public:
      /*! \brief Default constructor. */
      EvolutionParameters();
      
      /*! \brief Cloning operator. */
      EvolutionParameters<R>* clone() const;
      
      /*! \brief The maximum number of steps for an iterative algorithm. */
      size_type maximum_number_of_steps() const;

      /*! \brief The time after which an MapEvolver may approximate computed sets on a grid,
       *  in order to use previously cached results for the grid. Increasing this 
       *  parameter may improve the accuracy of the computations.  
       *  If there is recurrence in the system, then this parameter should be set to 
       *  the average recurrence time, if known. Used for discrete-time computation.
       */
      size_type lock_to_grid_steps() const;


      /*! \brief A suggested minimum step size for integration. This value may be ignored if an integration step cannot be performed without reducing the step size below this value. */
      time_type minimum_step_size() const;
      /*! \brief The maximum allowable step size for integration. Decreasing this value increases the accuracy of the computation. */
      time_type maximum_step_size() const;
      /*! \brief A suggested minimum radius of a basic set after a subdivision (not a strict bound). */
      R minimum_basic_set_radius() const;
      /*! \brief The maximum allowable radius of a basic set during integration. Decreasing this value increases the accuracy of the computation of an over-approximation. */
      R maximum_basic_set_radius() const;

      /*! \brief The time after which an VectorFieldEvolver may approximate computed sets on a grid,
       *  in order to use previously cached integration results for the grid. Increasing this 
       *  parameter improves the accuracy of the computations. Setting this parameter too
       *  low usually results in meaningless computations. A typical system trajectory 
       *  should move at least four times the grid size between locking to the grid. <br>
       *  For forced oscillators, this parameter should be set to the forcing time, 
       *  or a multiple or fraction thereof. Used for continuous-time computation.
       * 
       */
      time_type lock_to_grid_time() const;

      /*! \brief Set the length of the approximation grid. Decreasing this value increases the accuracy of the computation. */
      R grid_length() const;
      /*! \brief Set the default length of the approximation grid used for the argument of the function. Decreasing this value increases the accuracy of the computation.  */
      R argument_grid_length() const;
      /*! \brief Set the default length of the approximation grid for the result of the function. Decreasing this value increases the accuracy of the computation.  */
      R result_grid_length() const;

      /*! \brief Set the size of the region used for computation. Increasing this value reduces the risk of error due to missing orbits which leave the bounding domain. */
      R bounding_domain_size() const;

      /*! \brief A grid of dimension \a d with the default spacing. */
      Geometry::Rectangle<R> bounding_box(dimension_type d) const;

      /*! \brief A grid of dimension \a d with the default spacing. */
      Geometry::Grid<R> grid(dimension_type d) const;

      /*! \brief A grid of dimension \a d with the default spacing and bounds. */
      Geometry::FiniteGrid<R> finite_grid(dimension_type d) const;



      /*! \brief Set the maximum number of steps for an iterative algorithm. */
      void set_maximum_number_of_steps(size_type);
      /*! \brief Set the number of steps after which an applicator may approximate computed sets on a grid. */
      void set_lock_to_grid_steps(size_type);


      /*! \brief Set the suggested minimum step size for integration. */
      void set_minimum_step_size(time_type);
      /*! \brief Set the suggested maximum allowable step size for integration. */
      void set_maximum_step_size(time_type);
      
      /*! \brief Set the minimum radius of a basic set after a subdivision. */
      void set_minimum_basic_set_radius(R);
      /*! \brief Set the maximum radius of a basic set after a subdivision. */
      void set_maximum_basic_set_radius(R);

      /*! \brief Set the time after which an integrator may approximate computed sets on a grid. */
      void set_lock_to_grid_time(time_type);

      /*! \brief Set the length of the approximation grid. */
      void set_grid_length(R);
      /*! \brief Set the length of the approximation grid for the argument of a function. */
      void set_argument_grid_length(R);
      /*! \brief Set the length of the approximation grid for the result of a function. */
      void set_result_grid_length(R);

      /*! \brief Set the size of the region used for computation. */
      void set_bounding_domain_size(R);

      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;      
    };

    template<class R> inline 
    std::ostream& operator<<(std::ostream& os, const EvolutionParameters<R>& p) { 
      return p.write(os);
    }

  }
}

#endif /* ARIADNE_EVOLUTION_PARAMETERS_H */
