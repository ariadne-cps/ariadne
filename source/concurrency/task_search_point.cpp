/***************************************************************************
 *            concurrency/task_search_point.cpp
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

#include "task_search_point.hpp"
#include "task_search_space.hpp"
#include "../output/logging.hpp"

namespace Ariadne {

TaskSearchPoint::TaskSearchPoint(TaskSearchSpace const& space, ParameterBindingsMap const& bindings)
    : _space(space.clone()), _bindings(bindings) { }

TaskSearchPoint::TaskSearchPoint(TaskSearchPoint const& p) {
    this->_bindings.clear();
    this->_bindings.adjoin(p._bindings);
    this->_CACHED_SHIFT_BREADTHS = p._CACHED_SHIFT_BREADTHS;
    this->_space.reset(p.space().clone());
}

Set<TaskSearchPoint> TaskSearchPoint::make_random_shifted(Nat amount) const {
    Set<TaskSearchPoint> result;
    TaskSearchPoint current_point = *this;
    result.insert(current_point);
    while (result.size() < amount) {
        List<Nat> breadths = current_point.shift_breadths();
        Nat total_breadth = 0; for (auto b : breadths) total_breadth += b;
        auto space = this->space();

        Nat offset = (Nat)rand() % total_breadth;
        Nat current_breadth = 0;
        ParameterBindingsMap shifted_bindings;
        Bool shifted = false;
        for (auto binding : current_point.bindings()) {
            auto const& param = current_point.parameter(binding.first);
            int value = binding.second;
            current_breadth += breadths.at(space.index(param));
            if (not shifted and current_breadth > offset) {
                value = param.shifted_value_from(value);
                shifted = true;
            }
            shifted_bindings.insert(std::pair<Identifier,int>(param.name(), value));
        }
        result.insert(space.make_point(shifted_bindings));

        Nat new_choice = (Nat)rand() % result.size();
        auto iter = result.begin();
        for (Nat i=0; i<new_choice; ++i) ++iter;
        current_point = *iter;
    }
    return result;
}

TaskSearchPoint TaskSearchPoint::make_adjacent_shifted() const {
    List<Nat> breadths = this->shift_breadths();
    Nat total_breadth = 0;
    for (auto b : breadths) total_breadth += b;
    ARIADNE_PRECONDITION(total_breadth != 0);
    Set<TaskSearchPoint> result;
    auto space = this->space();
    Nat offset = (Nat)rand() % total_breadth;

    Nat current_breadth = 0;
    ParameterBindingsMap shifted_bindings;
    Bool shifted = false;
    for (auto binding : _bindings) {
        auto const &param = parameter(binding.first);
        int value = binding.second;
        current_breadth += breadths.at(space.index(param));
        if (not shifted and current_breadth > offset) {
            value = param.shifted_value_from(value);
            shifted = true;
        }
        shifted_bindings.insert(std::pair<Identifier, int>(param.name(), value));
    }
    return space.make_point(shifted_bindings);
}

TaskSearchSpace const& TaskSearchPoint::space() const {
    return *_space;
}

List<int> TaskSearchPoint::coordinates() const {
    return _bindings.values();
}

ParameterBindingsMap const& TaskSearchPoint::bindings() const {
    return _bindings;
}

int TaskSearchPoint::value(Identifier const& name) const {
    return _bindings.at(name);
}

SizeType TaskSearchPoint::index(Identifier const& name) const {
    return _space->index(name);
}

TaskSearchParameter const& TaskSearchPoint::parameter(Identifier const& name) const {
    return _space->parameter(name);
}

TaskSearchPoint& TaskSearchPoint::operator=(TaskSearchPoint const& p) {
    this->_bindings.clear();
    this->_bindings.adjoin(p._bindings);
    this->_CACHED_SHIFT_BREADTHS = p._CACHED_SHIFT_BREADTHS;
    this->_space.reset(p.space().clone());
    return *this;
}

Bool TaskSearchPoint::operator==(TaskSearchPoint const& p) const {
    for (auto iter=_bindings.begin(); iter!=_bindings.end(); ++iter) {
        if (p._bindings.at(iter->first) != iter->second)
            return false;
    }
    return true;
}

Bool TaskSearchPoint::operator<(TaskSearchPoint const& p) const {
    for (auto b : _bindings) {
        auto const this_value = b.second;
        auto const other_value = p._bindings.at(b.first);
        if (this_value < other_value) return true;
        else if (this_value > other_value) return false;
    }
    return false; // They are equal
}

Nat TaskSearchPoint::distance(TaskSearchPoint const& p) const {
    Nat result = 0;
    for (auto b : _bindings) {
        auto const& param = parameter(b.first);
        auto const v1 = b.second;
        auto const v2 = p._bindings.at(param.name());
        if (param.is_metric()) result += (v1 > v2 ? (Nat)(v1 - v2) : (Nat)(v2 - v1));
        else result += (v1 == v2 ? 0 : 1);
    }
    return result;
}

OutputStream& TaskSearchPoint::_write(OutputStream& os) const {
    return os << _bindings.values();
}

List<Nat> TaskSearchPoint::shift_breadths() const {
    if (_CACHED_SHIFT_BREADTHS.empty()) {
        for (auto b : _bindings) {
            auto const& param = parameter(b.first);
            auto const& values = param.values();
            auto size = values.size();
            if (not param.is_metric()) _CACHED_SHIFT_BREADTHS.append(size-1); // all except the current
            else if (b.second == values[size-1]) _CACHED_SHIFT_BREADTHS.append(1); // can only move down
            else if (b.second == values[0]) _CACHED_SHIFT_BREADTHS.append(1); // can only move up
            else _CACHED_SHIFT_BREADTHS.append(2);; // can move either up or down
        }
    }
    return _CACHED_SHIFT_BREADTHS;
}

Set<TaskSearchPoint> make_extended_set_by_shifting(Set<TaskSearchPoint> const& sources, SizeType size) {
    ARIADNE_PRECONDITION(size>=sources.size());
    ARIADNE_PRECONDITION(sources.begin()->space().total_points() >= size);
    auto expanded_sources = sources; // To be be expanded if the previous sources are incapable of getting the required size
    auto result = sources;

    while (result.size() < size) {
        auto source_it = expanded_sources.begin();
        SizeType previous_size = result.size();
        while (result.size() < size) {
            result.insert(source_it->make_adjacent_shifted());
            ++source_it; // Will move to next source even if no shift has been found
            if (source_it == expanded_sources.end()) break;
        }
        // At the end of a cycle, if no increment of points has been obtained, add the current results to the sources
        if (result.size() == previous_size) expanded_sources.adjoin(result);
    }
    return result;
}

} // namespace Ariadne
