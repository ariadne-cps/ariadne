/***************************************************************************
 *            concurrency/task_search_space.hpp
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

/*! \file concurrency/task_search_space.hpp
 *  \brief Class for handling a search space of parameters for a task.
 */

#ifndef ARIADNE_TASK_SEARCH_SPACE_HPP
#define ARIADNE_TASK_SEARCH_SPACE_HPP

#include "../utility/writable.hpp"
#include "../concurrency/task_search_parameter.hpp"

namespace Ariadne {

class TaskSearchPoint;
class Real;
template<class R> class Variable;

using ParameterBindingsMap = Map<TaskSearchParameter,int>;

class TaskSearchSpace : public WritableInterface {
  public:
    TaskSearchSpace(Set<TaskSearchParameter> const& parameters);

    TaskSearchPoint make_point(Map<Identifier,int> const& bindings) const;
    TaskSearchPoint make_point(ParameterBindingsMap const& bindings) const;
    TaskSearchPoint initial_point() const;

    List<TaskSearchParameter> const& parameters() const;

    //! \brief The total number of points identified by the space
    SizeType total_points() const;
    //! \brief The number of parameters in the space
    SizeType dimension() const;
    //! \brief If there are no parameters in the space
    Bool is_empty() const;
    //! \brief The index of the given parameter in the ordered space
    SizeType index(TaskSearchParameter const& p) const;

    TaskSearchSpace* clone() const;

    virtual OutputStream& _write(OutputStream& os) const;

  private:
    const List<TaskSearchParameter> _parameters;
};

} // namespace Ariadne

#endif // ARIADNE_TASK_SEARCH_SPACE_HPP
