/***************************************************************************
 *            dynamics/simulation_toolbox_interface.hpp
 *
 *  Copyright  2008-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file dynamics/simulation_toolbox_interface.hpp
 *  \brief Interfaces for calculus tools useful for working with dynamical systems.
 */


#ifndef ARIADNE_SIMULATION_TOOLBOX_INTERFACE_HPP
#define ARIADNE_SIMULATION_TOOLBOX_INTERFACE_HPP

#include "../utility/tribool.hpp"

namespace Ariadne {

class ScalarMultivariateFunction;
class VectorMultivariateFunction;

/*! \brief Tools for analysing dynamical systems by simulation.
 */
class SimulationToolboxInterface
{
  public:
    typedef FloatDP RealType;
    typedef FloatDP TimeType;
    typedef ExactPoint StateType;
    typedef ScalarMultivariateFunction PredicateType;
    typedef VectorMultivariateFunction MapType;
    typedef VectorMultivariateFunction VectorFieldType;
  public:
    //! \brief Virtual destructor.
    virtual ~SimulationToolboxInterface() = default;

    //! \brief Test if a box satisfies the constraint given by the guard. Returns \a true is the point
    //! satisfies the constraint, \a false if the point does not satisfy the constraint, and
    //! indeterminate if the point lies on the boundary.
    virtual ValidatedKleenean
    active(const ScalarMultivariateFunction& guard,
           const ExactPoint& state) const = 0;

    //! \brief Computes the time at which points in the \a initial_point cross the zero-set of the
    //! the \a guard under evolution of the \a vector_field, for times up to \a maximum_time.
    //! The crossing must be (differentiably) transverse.
    virtual TimeType
    crossing_time(const ScalarMultivariateFunction& guard,
                  const VectorMultivariateFunction& vector_field,
                  const ExactPoint& initial_state,
                  const TimeType& maximum_time) const = 0;

    //! \brief Computes the image of the set defined by \a state under the \a map.
    virtual StateType
    reset_step(const MapType& map,
               const ExactPoint& state) const = 0;

    //! \brief Computes the points reached by evolution of the \a initial_state under the flow
    //! given by \a vector_field. The \a step_size gives the time the points should be flowed.
    virtual StateType
    integration_step(const VectorMultivariateFunction& vector_field,
                     const ExactPoint& initial_state,
                     const TimeType& step_size) const = 0;

};

} //  namespace Ariadne


#endif // ARIADNE_SIMULATION_TOOLBOX_INTERFACE_HPP */
