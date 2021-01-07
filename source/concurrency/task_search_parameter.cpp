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

#include "task_search_parameter.hpp"

namespace Ariadne {

std::ostream& operator<<(std::ostream& os, const TaskSearchParameterKind kind) {
    switch(kind) {
    case TaskSearchParameterKind::BOOLEAN:
        os << "Boolean";
        break;
    case TaskSearchParameterKind::ENUMERATION:
        os << "Enumeration";
        break;
    case TaskSearchParameterKind::METRIC:
        os << "Metric";
        break;
    default:
        ARIADNE_FAIL_MSG("Unhandled value " << kind << " in TaskSearchParameterKind for printing.");
    }
    return os;
}

Real TaskSearchParameterBase::value(Nat integer_value, Map<RealVariable,Real> const& constant_values) const {
    Map<Identifier,Real> conversion_variables;
    for (auto entry : constant_values) conversion_variables[entry.first.name()] = entry.second;
    conversion_variables.insert(Pair<Identifier,Real>(_variable.name(),integer_value));

    return evaluate(_value_expression, Valuation<Real,Real>(conversion_variables));
}

Nat MetricSearchParameter::shifted_value_from(Nat value) const {
    if (value == 0) return 1;
    if (value == upper_bound()) return value-1;
    if (rand() % 2 == 0) return value-1;
    else return value+1;
}

Nat BooleanSearchParameter::shifted_value_from(Nat value) const {
    return (value == 1 ? 0 : 1);
}

Bool TaskSearchParameter::operator==(TaskSearchParameter const& p) const {
    // FIXME: enumeration tasks should be distinguished
    return _impl->name() == p._impl->name() and _impl->kind() == p._impl->kind() and _impl->upper_bound() == p._impl->upper_bound();
}

Bool TaskSearchParameter::operator<(TaskSearchParameter const& p) const {
    return name() < p.name();
}

OutputStream& TaskSearchParameter::_write(OutputStream& os) const {
    return os << "('" << name() << "', kind: " << kind() << ", upper_bound: " + to_string(upper_bound()) + ", initial: " << to_string(initial()) << ")";
}

} // namespace Ariadne
