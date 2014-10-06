/***************************************************************************
 *            simulation_toolbox_interface.h
 *
 *  Copyright  2008  Pieter Collins
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

/*! \file simulation_toolbox_interface.h
 *  \brief Interfaces for calculus tools useful for working with dynamical systems.
 */


#ifndef ARIADNE_SIMULATION_TOOLBOX_INTERFACE_H
#define ARIADNE_SIMULATION_TOOLBOX_INTERFACE_H

#include "tribool.h"

namespace Ariadne {

class ScalarFunction;
class VectorFunction;

/*! \brief Tools for analysing dynamical systems by simulation.
 */
class SimulationToolboxInterface
{
  public:
    typedef Float RealType;
    typedef Float TimeType;
    typedef ExactPoint StateType;
    typedef ScalarFunction PredicateType;
    typedef VectorFunction MapType;
    typedef VectorFunction VectorFieldType;
  public:
    //! \brief Virtual destructor.
    virtual ~SimulationToolboxInterface() { }

    //! \brief Test if a box satisfies the constraint given by the guard. Returns \a true is the point
    //! satisfies the constraint, \a false if the point does not satisfy the constraint, and
    //! indeterminate if the point lies on the boundary.
    virtual tribool
    active(const ScalarFunction& guard,
           const ExactPoint& state) const = 0;

    //! \brief Computes the time at which points in the \a initial_point cross the zero-set of the
    //! the \a guard under evolution of the \a vector_field, for times up to \a maximum_time.
    //! The crossing must be (differentiably) transverse.
    virtual TimeType
    crossing_time(const ScalarFunction& guard,
                  const VectorFunction& vector_field,
                  const ExactPoint& initial_state,
                  const TimeType& maximum_time) const = 0;

    //! \brief Computes the image of the set defined by \a state under the \a map.
    virtual StateType
    reset_step(const MapType& map,
               const ExactPoint& state) const = 0;

    //! \brief Computes the points reached by evolution of the \a initial_state under the flow
    //! given by \a vector_field. The \a step_size gives the time the points should be flowed.
    virtual StateType
    integration_step(const VectorFunction& vector_field,
                     const ExactPoint& initial_state,
                     const TimeType& step_size) const = 0;

};

} //  namespace Ariadne


#endif // ARIADNE_SIMULATION_TOOLBOX_INTERFACE_H */
