/***************************************************************************
 *            concurrency/task_appraisal.cpp
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

#include "../concurrency/task_appraisal.hpp"

namespace Ariadne {

TaskSearchPointAppraisal::TaskSearchPointAppraisal(TaskSearchPoint const& p,
                                                   CostType const& s,
                                                   SizeType const& permissive_failures,
                                                   SizeType const& critical_failures)
           : _point(p), _cost(s), _permissive_failures(permissive_failures), _critical_failures(critical_failures) { }

TaskSearchPoint const&
TaskSearchPointAppraisal::point() const {
    return _point;
}

CostType const&
TaskSearchPointAppraisal::cost() const {
    return _cost;
}

SizeType const&
TaskSearchPointAppraisal::permissive_failures() const {
    return _permissive_failures;
}

SizeType const&
TaskSearchPointAppraisal::critical_failures() const {
    return _critical_failures;
}

Bool TaskSearchPointAppraisal::operator<(TaskSearchPointAppraisal const& s) const {
    if (this->_critical_failures > s._critical_failures) return false;
    else if (this->_critical_failures < s._critical_failures) return true;
    else if (this->_permissive_failures > s._permissive_failures) return false;
    else if (this->_permissive_failures < s._permissive_failures) return true;
    else return this->_cost < s._cost;
}

OutputStream& TaskSearchPointAppraisal::_write(OutputStream& os) const {
    return os << "{" << _point << ":" << _cost << (_permissive_failures > 0 ? (",P:" + to_string(_permissive_failures)) : "")
              << (_critical_failures > 0 ? (",C:" + to_string(_critical_failures)) : "") << "}";
}
} // namespace Ariadne
