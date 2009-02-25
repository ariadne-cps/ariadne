/***************************************************************************
 *            simulation_toolbox.h
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
 
/*! \file simulation_toolbox.h
 *  \brief Methods of point calculus based on the PointVariable class.
 */


#ifndef ARIADNE_SIMULATION_TOOLBOX_H
#define ARIADNE_SIMULATION_TOOLBOX_H

#include "tribool.h"
#include "logging.h"
#include "simulation_toolbox_interface.h"

/* \brief Top-level namespace. */
namespace Ariadne {
 

typedef double Float;
class ExpressionInterface;
class FunctionInterface;
class Point;

typedef std::pair<Point,Point> Segment;

/*! \brief Tools for analysing dynamical systems based on approximate simulation. */
class SimulationToolbox
    : public SimulationToolboxInterface
{
  public:
    typedef Float RealType;
    typedef Float TimeType;
  public:
    //! \brief Default constructor.
    SimulationToolbox();

    //! \brief Test if a set satisfied the constraint given by the guard model. Returns \a true is all 
    //! points in the set satisfy the constraint, \a false if all points do not satisfy the constraint, 
    //! and indeterminate otherwise.
    tribool active(const ExpressionInterface& guard, 
                   const Point& point) const;

  
    //! \brief Computes the time at which point  \a initial_point cross the zero-set of the
    //! the \a guard under evolution of the \a vector_field, for times up to \a maximum_time.
    //! The crossing must be (differentiably) transverse.
    TimeType crossing_time(const ExpressionInterface& guard,
                           const FunctionInterface& vector_field, 
                           const Point& initial_point, 
                           const TimeType& maximum_time) const;

    //! \brief Computes the image of the set defined by \a point under the map \a map. 
    Point reset_step(const FunctionInterface& map, 
                     const Point& point) const;
  
    //! \brief Computes the points reached by evolution of the \a initial_point under the flow
    //! given by the \a vector_field. The \a step_size gives the time the point 
    //! should be flowed.
    Point integration_step(const FunctionInterface& vector_field, 
                           const Point& initial_point, 
                           const TimeType& step_size) const;
  

};

}


#endif /* ARIADNE_SIMULATION_TOOLBOX_H */
