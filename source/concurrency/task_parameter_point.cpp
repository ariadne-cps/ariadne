/***************************************************************************
 *            concurrency/task_parameter_point.cpp
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

TaskParameterPoint::TaskParameterPoint(TaskParameterSpace const& space, ParameterBindingsMap const& bindings)
    : _space(space.clone()), _bindings(bindings) { }

Set<TaskParameterPoint> TaskParameterPoint::make_random_shifted(Nat amount) const {
    Set<TaskParameterPoint> result;
    TaskParameterPoint current_point = *this;
    for (Nat num_points=1; num_points<=amount; ++num_points) {
        List<Nat> breadths = current_point.shift_breadths();
        Nat total_breadth = 0;
        for (Nat b : breadths) total_breadth += b;
        auto space = this->space();

        while(true) {
            Nat offset = (Nat)rand() % total_breadth;
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
    do offsets.insert((Nat)rand() % total_breadth); while (offsets.size() < amount);

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

Map<TaskParameter,Real> TaskParameterPoint::values(Map<RealVariable,Real> const & external_values) const {
    Map<TaskParameter,Real> result;
    for (auto binding : _bindings) {
        result.insert(Pair<TaskParameter,Real>(binding.first,binding.first.value(binding.second,external_values)));
    }
    return result;
}

Real TaskParameterPoint::value(Identifier const& var, Map<RealVariable,Real> const& external_values) const {
    auto iter = _bindings.begin();
    for (;iter!=_bindings.end();++iter) {
        if (iter->first.name() == var) break;
    }
    if (iter == _bindings.end())
        ARIADNE_FAIL_MSG("Variable '" << var << "' could not be found in the space.");
    return iter->first.value(iter->second,external_values);
}

Real TaskParameterPoint::time_cost_estimate(Map<RealVariable,Real> const& external_variables) const {
    Map<Identifier,Real> values;
    for (auto binding : _bindings) {
        values.insert(Pair<Identifier,Real>(binding.first.variable().name(),binding.first.value(binding.second,external_variables)));
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
        if (iter->first.kind() == TaskParameterKind::METRIC) result += (v1 > v2 ? v1 - v2 : v2 - v1);
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
            if (iter->first.kind() != TaskParameterKind::METRIC) _CACHED_SHIFT_BREADTHS.append(iter->first.upper_bound()); // size-1 choices
            else if (iter->second == iter->first.upper_bound()) _CACHED_SHIFT_BREADTHS.append(1); // can only move down
            else if (iter->second == 0) _CACHED_SHIFT_BREADTHS.append(1); // can only move up
            else _CACHED_SHIFT_BREADTHS.append(2);; // can move either up or down
        }
    }
    return _CACHED_SHIFT_BREADTHS;
}

Set<TaskParameterPoint> make_adjacent_set_shifted_from(Set<TaskParameterPoint> const& sources, Nat amount) {
    Set<TaskParameterPoint> result;
    result.adjoin(sources);
    for (TaskParameterPoint source: sources) {
        SizeType previous_size = result.size();
        Nat remaining_amount = amount;
        do {
            result.adjoin(source.make_adjacent_shifted(amount));
            remaining_amount = amount-(result.size()-previous_size);
            previous_size = result.size();
        } while (remaining_amount > 0);
    }
    return result;
}

} // namespace Ariadne
