/***************************************************************************
 *            dynamics/task_parameter.cpp
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

#include "dynamics/task_parameter.hpp"

namespace Ariadne {

bool unique_elements(const std::vector<TaskParameter>& lst) {
    Set<TaskParameter> found;
    for(auto iter=lst.begin(); iter!=lst.end(); ++iter) {
        if(found.contains(*iter)) { return false; } else { found.insert(*iter); } }
    return true;
}

template <typename E>
constexpr typename std::underlying_type<E>::type to_underlying(E e) {
    return static_cast<typename std::underlying_type<E>::type>(e);
}

Nat MetricTaskParameter::shifted_value_from(Nat value) const {
    if (value == 0) return 1;
    if (is_upper_bounded() and value == upper_bound()) return value-1;
    if (rand() % 2 == 0) return value-1;
    else return value+1;
}

Nat BooleanTaskParameter::shifted_value_from(Nat value) const {
    return (value == 1 ? 0 : 1);
}

Bool TaskParameter::operator==(TaskParameter const& p) const {
    return _impl->name() == p._impl->name() and _impl->is_metric() == p._impl->is_metric() and _impl->upper_bound() == p._impl->upper_bound();
}

Bool TaskParameter::operator<(TaskParameter const& p) const {
    return name() < p.name();
}

OutputStream& TaskParameter::_write(OutputStream& os) const {
    return os << "(" << name() << ", " << is_metric() << (is_upper_bounded() ? (", " + to_string(upper_bound()) + ")") : ")");
}

TaskParameterSpace::TaskParameterSpace(List<TaskParameter> const& parameters) : _parameters(parameters) {
    ARIADNE_PRECONDITION(not _parameters.empty());
    ARIADNE_PRECONDITION(unique_elements(parameters));
}

Nat TaskParameterSpace::index(TaskParameter const& p) const {
    for (SizeType i=0; i<_parameters.size(); ++i) if (_parameters.at(i) == p) return i;
    ARIADNE_FAIL_MSG("Task parameter '" << p << "' not found in the space.");
}

OutputStream& TaskParameterSpace::_write(OutputStream& os) const {
    return os << _parameters;
}


TaskParameterSpace TaskParameterPoint::space() const {
    return TaskParameterSpace(make_list(_bindings.keys()));
}

//! \brief Generate a new point by shifting a given \a amount
//! \details The amount is the total number of shifts allowed across the parameters
TaskParameterPoint TaskParameterPoint::shift(Nat amount) const {
    List<Nat> breadths = _shift_breadths();
    Nat total_breadth = 0;
    for (Nat b : breadths) total_breadth += b;
    Set<Nat> offsets;
    do offsets.insert((Nat)rand() % total_breadth); while (offsets.size() < amount);

    ParameterBindingsMap shifted_bindings;
    Nat current_breadth = 0;
    auto breadth_iter = breadths.begin();
    auto offset_iter = offsets.begin();
    for (auto param_iter=_bindings.begin(); param_iter != _bindings.end(); ++param_iter) {
        Nat val = param_iter->second;
        if (offset_iter != offsets.end() and current_breadth >= *offset_iter) { // If this should be shifted
            val = param_iter->first.shifted_value_from(param_iter->second); // Shift
            ++offset_iter;
        }
        current_breadth += *breadth_iter;
        ++breadth_iter;
        shifted_bindings.insert(std::pair<TaskParameter,Nat>(param_iter->first, val));
    }
    return TaskParameterPoint(shifted_bindings);
}

Bool TaskParameterPoint::operator==(TaskParameterPoint const& p) const {
    for (auto iter=_bindings.begin(); iter!=_bindings.end(); ++iter) {
        if (p._bindings.at(iter->first) != iter->second)
            return false;
    }
    return true;
}

OutputStream& TaskParameterPoint::_write(OutputStream& os) const {
    return os << _bindings.values();
}

List<Nat> TaskParameterPoint::_shift_breadths() const {
    if (_CACHED_SHIFT_BREADTHS.empty()) {
        for (auto iter = _bindings.begin(); iter != _bindings.end(); ++iter) {
            if (not iter->first.is_metric()) _CACHED_SHIFT_BREADTHS.append(iter->first.upper_bound()); // size-1 choices
            else if (iter->first.is_upper_bounded() and iter->second == iter->first.upper_bound()) _CACHED_SHIFT_BREADTHS.append(1); // can only move down
            else if (iter->second == 0) _CACHED_SHIFT_BREADTHS.append(1); // can only move up
            else _CACHED_SHIFT_BREADTHS.append(2);; // can move either up or down
        }
    }
    return _CACHED_SHIFT_BREADTHS;
}

} // namespace Ariadne
