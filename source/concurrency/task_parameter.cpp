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

#include "task_parameter.hpp"

namespace Ariadne {

template <typename E>
constexpr typename std::underlying_type<E>::type to_underlying(E e) {
    return static_cast<typename std::underlying_type<E>::type>(e);
}

Nat MetricTaskParameter::shifted_value_from(Nat value) const {
    if (value == 0) return 1;
    if (value == upper_bound()) return value-1;
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
    return os << "(" << name() << ", metric? " << is_metric() << ", ub: " + to_string(upper_bound()) + ")";
}

TaskParameterSpace::TaskParameterSpace(Set<TaskParameter> const& parameters, RealExpression const& time_cost_estimator)
: _parameters(parameters), _time_cost_estimator(time_cost_estimator) {
    ARIADNE_PRECONDITION(not _parameters.empty());
}

TaskParameterPoint TaskParameterSpace::make_point(ParameterBindingsMap const& bindings) const {
    return TaskParameterPoint(*this,bindings);
}

TaskParameterPoint TaskParameterSpace::make_point(Map<RealVariable,Nat> const& bindings) const {
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

Set<TaskParameterPoint> TaskParameterPoint::make_random_shifted(Nat amount) const {
    Set<TaskParameterPoint> result;
    TaskParameterPoint current_point = *this;
    for (Nat num_points=1; num_points<=amount; ++num_points) {
        List<Nat> breadths = current_point.shift_breadths();
        Nat total_breadth = 0;
        for (Nat b : breadths) total_breadth += b;
        auto space = this->space();

        while(true) {
            Nat offset = (Nat)rand() % (total_breadth-1);
            Nat current_breadth = 0;
            ParameterBindingsMap shifted_bindings;
            Bool shifted = false;
            for (auto binding : current_point.bindings()) {
                auto const& param = binding.first;
                Nat value = binding.second;
                current_breadth += breadths.at(space.index(param));
                if (not shifted and current_breadth > offset) {
                    value = param.shifted_value_from(binding.second);
                    shifted = true;
                }
                shifted_bindings.insert(std::pair<TaskParameter,Nat>(param, value));
            }
            result.insert(space.make_point(shifted_bindings));

            Nat new_choice = (Nat)rand() % result.size();
            auto iter = result.begin();
            for (Nat i=1; i<new_choice; ++i) ++iter;
            current_point = *iter;

            if (result.size() == num_points) break;
        }
    }
    return result;
}

Set<TaskParameterPoint> TaskParameterPoint::make_adjacent_shifted(Nat amount) const {
    List<Nat> breadths = this->shift_breadths();
    Nat total_breadth = 0;
    for (Nat b : breadths) total_breadth += b;
    ARIADNE_PRECONDITION(total_breadth >= amount);
    Set<TaskParameterPoint> result;
    auto space = this->space();
    Set<Nat> offsets;
    do offsets.insert((Nat)rand() % (total_breadth-1)); while (offsets.size() < amount);

    Nat num_points = 1;
    for (Nat offset : offsets) {
        do {
            Nat current_breadth = 0;
            ParameterBindingsMap shifted_bindings;
            Bool shifted = false;
            for (auto binding : _bindings) {
                auto const& param = binding.first;
                Nat value = binding.second;
                current_breadth += breadths.at(space.index(param));
                if (not shifted and current_breadth > offset) {
                    value = param.shifted_value_from(binding.second);
                    shifted = true;
                }
                shifted_bindings.insert(std::pair<TaskParameter,Nat>(param, value));
            }
            result.insert(space.make_point(shifted_bindings));
        } while (result.size() < num_points);
        ++num_points;
    }
    return result;
}

Nat TaskParameterPoint::hash_code() const {
    Nat result=0;
    Nat prod=1;
    for (auto iter=_bindings.begin(); iter!=_bindings.end(); ++iter) {
        result += iter->second*prod;
        prod *= iter->first.upper_bound()+1;
    }
    return result;
}

Real TaskParameterPoint::time_cost_estimate() const {
    Map<Identifier,Real> values;
    for (auto binding : _bindings) {
        Map<Identifier,Real> conversion_variables;
        conversion_variables.insert(Pair<Identifier,Real>(binding.first.variable().name(),binding.second));
        Real converted_value = evaluate(binding.first.integer_conversion(), Valuation<Real,Real>(conversion_variables));
        values.insert(Pair<Identifier,Real>(binding.first.variable().name(),converted_value));
    }
    return evaluate(_space->time_cost_estimator(),Valuation<Real,Real>(values));
}

TaskParameterPoint& TaskParameterPoint::operator=(TaskParameterPoint const& p) {
    this->_bindings.clear();
    this->_bindings.adjoin(p._bindings);
    this->_CACHED_SHIFT_BREADTHS = p._CACHED_SHIFT_BREADTHS;
    this->_space.reset(p.space().clone());
    return *this;
}

Bool TaskParameterPoint::operator==(TaskParameterPoint const& p) const {
    for (auto iter=_bindings.begin(); iter!=_bindings.end(); ++iter) {
        if (p._bindings.at(iter->first) != iter->second)
            return false;
    }
    return true;
}

Bool TaskParameterPoint::operator<(TaskParameterPoint const& p) const {
    // ASSUMPTION: they have the same space
    for (auto iter=_bindings.begin(); iter!=_bindings.end(); ++iter) {
        Nat const this_value = iter->second;
        Nat const other_value = p._bindings.at(iter->first);
        if (this_value < other_value) return true;
        else if (this_value > other_value) return false;
    }
    return false; // They are equal
}

Nat TaskParameterPoint::distance(TaskParameterPoint const& p) const {
    // ASSUMPTION: they have the same space
    Nat result = 0;
    for (auto iter=_bindings.begin(); iter!=_bindings.end(); ++iter) {
        Nat const v1 = iter->second;
        Nat const v2 = p._bindings.at(iter->first);
        if (iter->first.is_metric()) result += (v1 > v2 ? v1 - v2 : v2 - v1);
        else result += (v1 == v2 ? 0 : 1);
    }
    return result;

}


OutputStream& TaskParameterPoint::_write(OutputStream& os) const {
    return os << _bindings.values();
}

List<Nat> TaskParameterPoint::shift_breadths() const {
    if (_CACHED_SHIFT_BREADTHS.empty()) {
        for (auto iter = _bindings.begin(); iter != _bindings.end(); ++iter) {
            if (not iter->first.is_metric()) _CACHED_SHIFT_BREADTHS.append(iter->first.upper_bound()); // size-1 choices
            else if (iter->second == iter->first.upper_bound()) _CACHED_SHIFT_BREADTHS.append(1); // can only move down
            else if (iter->second == 0) _CACHED_SHIFT_BREADTHS.append(1); // can only move up
            else _CACHED_SHIFT_BREADTHS.append(2);; // can move either up or down
        }
    }
    return _CACHED_SHIFT_BREADTHS;
}

} // namespace Ariadne
