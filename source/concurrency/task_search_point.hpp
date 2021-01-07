/***************************************************************************
 *            concurrency/task_search_point.hpp
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

/*! \file concurrency/task_search_point.hpp
 *  \brief Class for handling search points of tool parameters for a task.
 */

#ifndef ARIADNE_SEARCH_POINT_HPP
#define ARIADNE_SEARCH_POINT_HPP

#include "task_search_parameter.hpp"

namespace Ariadne {

class TaskSearchSpace;
using ParameterBindingsMap = Map<TaskSearchParameter,Nat>;

class TaskSearchPoint : public WritableInterface {
    friend class TaskSearchSpace;
  protected:
    TaskSearchPoint(TaskSearchSpace const& space, ParameterBindingsMap const& bindings);
  public:
    TaskSearchPoint(TaskSearchPoint const& p);
    ~TaskSearchPoint() = default;

    //! \brief The parameter space
    TaskSearchSpace const& space() const { return *_space; }

    //! \brief The coordinates in the natural space, according to the space ordering
    List<Nat> coordinates() const { return _bindings.values(); }
    //! \brief The values in the real space, possibly using some \a external_values
    Map<TaskSearchParameter,Real> values(Map<RealVariable,Real> const& external_values = Map<RealVariable,Real>()) const;
    //! \brief The value for a given parameter, possibly using some \a external values
    Real value(Identifier const& var, Map<RealVariable,Real> const& external_values = Map<RealVariable,Real>()) const;

    //! \brief Provide a time cost estimate, where the \a external_values bind constants to their values
    Real time_cost_estimate(Map<RealVariable,Real> const & external_values = Map<RealVariable,Real>()) const;

    //! \brief Generate an \a amount of new points by shifting one parameter each
    //! \details Guarantees that all points are different and with distance equal to 1
    Set<TaskSearchPoint> make_adjacent_shifted(Nat amount) const;
    //! \brief Generate an \a amount of new points by shifting one parameter each from the current point,
    //! then the next point to shift from is a random one from those already generated
    //! \details Guarantees that all points are different
    Set<TaskSearchPoint> make_random_shifted(Nat amount) const;

    ParameterBindingsMap const& bindings() const { return _bindings; }
    Nat const& operator[](TaskSearchParameter const& p) const { return _bindings.at(p); }

    TaskSearchPoint& operator=(TaskSearchPoint const& p);
    //! \brief Equality check is performed under the assumption that we always work with the same parameters,
    //! hence no space check is performed.
    Bool operator==(TaskSearchPoint const& p) const;
    //! \brief Ordering is based on point value
    Bool operator<(TaskSearchPoint const& p) const;
    //! \brief The distance with respect to another point
    //! \details Distance between values for non-metric parameters is either 1 or 0
    Nat distance(TaskSearchPoint const& p) const;

    //! \brief Compute the breadth of possible shifts of the point for each parameter
    List<Nat> shift_breadths() const;

    //! \brief A positive number uniquely identifying the point in the space
    Nat hash_code() const;

    virtual OutputStream& _write(OutputStream& os) const;
  private:

    SharedPointer<TaskSearchSpace> _space;
    ParameterBindingsMap _bindings;

    mutable List<Nat> _CACHED_SHIFT_BREADTHS;
};

//! \brief Generate an \a amount of new points from each point in \a sources, by shifting one parameter each
//! \return The original points plus the shifted ones
//! \details Guarantees that all points are different and with distance equal to 1 to the related source point
Set<TaskSearchPoint> make_adjacent_set_shifted_from(Set<TaskSearchPoint> const& sources, Nat amount);

} // namespace Ariadne

#endif // ARIADNE_SEARCH_POINT_HPP
