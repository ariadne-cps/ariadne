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

#include "concurrency/task_search_parameter.hpp"

namespace Ariadne {

TaskSearchParameter::TaskSearchParameter(Identifier const& name, Bool is_metric, List<int> const& values) :
    _name(name), _is_metric(is_metric), _values(values) {
    ARIADNE_PRECONDITION(values.size()>1);
}

Identifier const& TaskSearchParameter::name() const {
    return _name;
}

List<int> const& TaskSearchParameter::values() const {
    return _values;
}

Bool TaskSearchParameter::is_metric() const {
    return _is_metric;
}

int TaskSearchParameter::random_value() const {
    return _values[(SizeType)rand() % _values.size()];
}

int TaskSearchParameter::shifted_value_from(int value) const {
    SizeType num_values = _values.size();
    if (_is_metric) {
        if (value == _values[0]) return value+1;
        if (value == _values[num_values-1]) return value-1;
        if ((SizeType)rand() % 2 == 0) return value+1;
        else return value-1;
    } else {
        SizeType rand_value = (SizeType)rand() % (num_values-1);
        if (_values[rand_value] == value) return _values[num_values-1];
        else return _values[rand_value];
    }
}

Bool TaskSearchParameter::operator==(TaskSearchParameter const& p) const {
    return name() == p.name();
}

Bool TaskSearchParameter::operator<(TaskSearchParameter const& p) const {
    return name() < p.name();
}

OutputStream& TaskSearchParameter::_write(OutputStream& os) const {
    os << "{'" << name() << "', is_metric=" << _is_metric << ", values=";
    if (_is_metric) os << "[" << _values[0] << ":" << _values[_values.size()-1] << "]";
    else os << _values;
    return os << "}";
}

} // namespace Ariadne
