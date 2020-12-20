/***************************************************************************
 *            concurrency/task_parameter_space.cpp
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

#include "task_parameter_point.hpp"
#include "task_parameter_space.hpp"

namespace Ariadne {

TaskParameterSpace::TaskParameterSpace(Set<TaskParameter> const& parameters, RealExpression const& time_cost_estimator)
        : _parameters(parameters), _time_cost_estimator(time_cost_estimator) {
    ARIADNE_PRECONDITION(not _parameters.empty());
}

TaskParameterPoint TaskParameterSpace::make_point(ParameterBindingsMap const& bindings) const {
    ARIADNE_PRECONDITION(bindings.size() == this->dimension())
    return TaskParameterPoint(*this,bindings);
}

TaskParameterPoint TaskParameterSpace::make_point(Map<RealVariable,Nat> const& bindings) const {
    ARIADNE_PRECONDITION(bindings.size() == this->dimension())
    ParameterBindingsMap pb;
    for (auto p : _parameters) {
        Nat v = bindings.find(p.variable())->second;
        pb.insert(Pair<TaskParameter,Nat>(p,v));
    }
    return TaskParameterPoint(*this,pb);
}

Nat TaskParameterSpace::index(TaskParameter const& p) const {
    for (SizeType i=0; i<_parameters.size(); ++i) if (_parameters.at(i) == p) return i;
    ARIADNE_FAIL_MSG("Task parameter '" << p << "' not found in the space.");
}

OutputStream& TaskParameterSpace::_write(OutputStream& os) const {
    return os << "{variables: "<< _parameters << ", time_cost_estimate: " << _time_cost_estimator << "}";
}

} // namespace Ariadne
