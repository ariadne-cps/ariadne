/***************************************************************************
 *            concurrency/task_parameter.cpp
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

} // namespace Ariadne
