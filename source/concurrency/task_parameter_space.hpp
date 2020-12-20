/***************************************************************************
 *            concurrency/task_parameter_space.hpp
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

/*! \file concurrency/task_parameter_space.hpp
 *  \brief Class for handling a space of tool parameters for a task.
 */

#ifndef ARIADNE_TASK_PARAMETER_SPACE_HPP
#define ARIADNE_TASK_PARAMETER_SPACE_HPP

namespace Ariadne {

class TaskParameter;
class TaskParameterPoint;

using ParameterBindingsMap = Map<TaskParameter,Nat>;

class TaskParameterSpace : public WritableInterface {
  public:
    TaskParameterSpace(Set<TaskParameter> const& parameters, RealExpression const& time_cost_estimator);

    TaskParameterPoint make_point(Map<RealVariable,Nat> const& bindings) const;
    TaskParameterPoint make_point(ParameterBindingsMap const& bindings) const;

    List<TaskParameter> const& parameters() const { return _parameters; }
    RealExpression const& time_cost_estimator() const { return _time_cost_estimator; }

    Nat dimension() const { return _parameters.size(); }
    Nat index(TaskParameter const& p) const;

    TaskParameterSpace* clone() const { return new TaskParameterSpace(*this); }

    virtual OutputStream& _write(OutputStream& os) const;

  private:
    const List<TaskParameter> _parameters;
    const RealExpression _time_cost_estimator;
};

} // namespace Ariadne

#endif // ARIADNE_TASK_PARAMETER_SPACE_HPP
