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
using ParameterBindingsMap = Map<ConfigurationPropertyPath,int>;

class TaskSearchPoint : public WritableInterface {
    friend class TaskSearchSpace;
  protected:
    TaskSearchPoint(TaskSearchSpace const& space, ParameterBindingsMap const& bindings);
  public:
    TaskSearchPoint(TaskSearchPoint const& p);
    ~TaskSearchPoint() = default;

    //! \brief The parameter space
    TaskSearchSpace const& space() const;

    //! \brief The coordinates in the natural space, according to the space ordering
    List<int> coordinates() const;

    //! \brief Generate a point adjacent to this one by shifting one parameter
    TaskSearchPoint make_adjacent_shifted() const;
    //! \brief Generate an \a amount of points by shifting one parameter each from the current point,
    //! then the next point to shift from is a random one from those already generated
    //! \details Guarantees that all points are different. Includes the original point.
    //! If \a amount is 1, no new point is generated.
    Set<TaskSearchPoint> make_random_shifted(Nat amount) const;

    ParameterBindingsMap const& bindings() const;

    //! \brief The value of the point for the given parameter path \a path
    int value(ConfigurationPropertyPath const& path) const;
    //! \brief The index in the space for the parameter identifier \a path
    SizeType index(ConfigurationPropertyPath const& path) const;
    //! \brief The parameter corresponding to the identifier \a path
    TaskSearchParameter const& parameter(ConfigurationPropertyPath const& path) const;

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

    virtual OutputStream& _write(OutputStream& os) const;
  private:

    SharedPointer<TaskSearchSpace> _space;
    ParameterBindingsMap _bindings;

    mutable List<Nat> _CACHED_SHIFT_BREADTHS;
};

//! \brief Generate an \a amount of new points from \a sources, by shifting one parameter each (ideally, see details)
//! \return The original points plus the shifted ones
//! \details \a size must be greater or equal than \a sources size but still lower than the maximum number of points
//! for the space. Shifting points are chosen by rotation, skipping to the next if the generated point is not new-
//! An effort is made to shift only by 1 with respect to the sources, but if not possible then the generated
//! points are added to the points used for shifting.
Set<TaskSearchPoint> make_extended_set_by_shifting(Set<TaskSearchPoint> const& sources, SizeType size);

} // namespace Ariadne

#endif // ARIADNE_SEARCH_POINT_HPP
