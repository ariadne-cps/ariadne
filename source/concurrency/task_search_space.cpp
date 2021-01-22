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

#include "concurrency/configuration_property_path.hpp"
#include "concurrency/task_search_point.hpp"
#include "concurrency/task_search_space.hpp"

namespace Ariadne {

TaskSearchSpace::TaskSearchSpace(Set<TaskSearchParameter> const& parameters)
        : _parameters(parameters) {
    ARIADNE_PRECONDITION(not _parameters.empty());
}

TaskSearchPoint TaskSearchSpace::make_point(ParameterBindingsMap const& bindings) const {
    ARIADNE_PRECONDITION(bindings.size() == this->dimension())
    ParameterBindingsMap pb;
    for (auto p : _parameters) {
        int v = bindings.find(p.path())->second;
        pb.insert(Pair<ConfigurationPropertyPath,int>(p.path(), v));
    }
    return TaskSearchPoint(*this, pb);
}

TaskSearchPoint TaskSearchSpace::initial_point() const {
    ParameterBindingsMap pb;
    for (auto p : _parameters) {
        pb.insert(Pair<ConfigurationPropertyPath,int>(p.path(), p.random_value()));
    }
    return TaskSearchPoint(*this, pb);
}

SizeType TaskSearchSpace::index(TaskSearchParameter const& p) const {
    for (SizeType i=0; i<_parameters.size(); ++i) if (_parameters.at(i) == p) return i;
    ARIADNE_FAIL_MSG("Task parameter '" << p << "' not found in the space.");
}

SizeType TaskSearchSpace::index(ConfigurationPropertyPath const& path) const {
    for (SizeType i=0; i<_parameters.size(); ++i) if (_parameters.at(i).path() == path) return i;
    ARIADNE_FAIL_MSG("Task parameter with path '" << path << "' not found in the space.");
}

TaskSearchParameter const& TaskSearchSpace::parameter(ConfigurationPropertyPath const& path) const {
    for (SizeType i=0; i<_parameters.size(); ++i) if (_parameters.at(i).path() == path) return _parameters.at(i);
    ARIADNE_FAIL_MSG("Task parameter with path '" << path << "' not found in the space.");
}

List<TaskSearchParameter> const& TaskSearchSpace::parameters() const {
    return _parameters;
}

SizeType TaskSearchSpace::total_points() const {
    SizeType result = 1;
    for (auto p : _parameters) result *= p.values().size();
    return result;
}

SizeType TaskSearchSpace::dimension() const {
    return _parameters.size();
}

TaskSearchSpace* TaskSearchSpace::clone() const {
    return new TaskSearchSpace(*this);
}

OutputStream& TaskSearchSpace::_write(OutputStream& os) const {
    return os << _parameters;
}

} // namespace Ariadne
