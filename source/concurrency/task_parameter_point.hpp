/***************************************************************************
 *            concurrency/task_parameter_point.hpp
 *
 *  Copyright  2007-20  Luca Geretti
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

/*! \file concurrency/task_parameter_point.hpp
 *  \brief Class for handling points of tool parameters for a task.
 */

#ifndef ARIADNE_TASK_PARAMETER_POINT_HPP
#define ARIADNE_TASK_PARAMETER_POINT_HPP

#include "task_parameter.hpp"

namespace Ariadne {

class TaskParameterSpace;
using ParameterBindingsMap = Map<TaskParameter,Nat>;

class TaskParameterPoint : public WritableInterface {
    friend class TaskParameterSpace;
  protected:
    TaskParameterPoint(TaskParameterSpace const& space, ParameterBindingsMap const& bindings);
  public:
    //! \brief The parameter space
    TaskParameterSpace const& space() const { return *_space; }

    //! \brief The values in the space, according to the space ordering
    List<Nat> values() const { return _bindings.values(); }

    //! \brief Provide a time cost estimate given the space's estimator
    Real time_cost_estimate() const;

    //! \brief Generate an \a amount of new points by shifting one parameter each
    //! \details Guarantees that all points are different and with distance equal to 1
    Set<TaskParameterPoint> make_adjacent_shifted(Nat amount) const;
    //! \brief Generate an \a amount of new points by shifting one parameter each from the current point,
    //! then the next point to shift from is a random one from those already generated
    //! \details Guarantees that all points are different
    Set<TaskParameterPoint> make_random_shifted(Nat amount) const;

    ParameterBindingsMap const& bindings() const { return _bindings; }
    Nat const& operator[](TaskParameter const& p) const { return _bindings.at(p); }

    TaskParameterPoint& operator=(TaskParameterPoint const& p);
    //! \brief Equality check is performed under the assumption that we always work with the same parameters,
    //! hence no space check is performed.
    Bool operator==(TaskParameterPoint const& p) const;
    //! \brief Ordering is based on point value
    Bool operator<(TaskParameterPoint const& p) const;
    //! \brief The distance with respect to another point
    //! \details Distance between values for non-metric parameters is either 1 or 0
    Nat distance(TaskParameterPoint const& p) const;

    //! \brief Compute the breadth of possible
    //! shifts of the point for each parameter
    List<Nat> shift_breadths() const;

    //! \brief A positive number uniquely identifying the point in the space
    Nat hash_code() const;

    virtual OutputStream& _write(OutputStream& os) const;
  private:

    SharedPointer<TaskParameterSpace> _space;
    ParameterBindingsMap _bindings;

    mutable List<Nat> _CACHED_SHIFT_BREADTHS;
};

//! \brief Generate an \a amount of new points from each point in \a sources, by shifting one parameter each
//! \return The original points plus the shifted ones
//! \details Guarantees that all points are different and with distance equal to 1 to the related source point
Set<TaskParameterPoint> make_adjacent_set_shifted_from(Set<TaskParameterPoint> const& sources, Nat amount);

} // namespace Ariadne

#endif // ARIADNE_TASK_PARAMETER_POINT_HPP
