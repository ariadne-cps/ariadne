/***************************************************************************
 *            concurrency/workload_advancement.cpp
 *
 *  Copyright  2007-21  Luca Geretti
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

#include "workload_advancement.hpp"

namespace Ariadne {

WorkloadAdvancement::WorkloadAdvancement(SizeType initial) : _num_waiting(initial), _num_processing(0), _num_completed(0) { }

SizeType WorkloadAdvancement::waiting() const {
    LockGuard<Mutex> lock(_mux);
    return _num_waiting;
}

SizeType WorkloadAdvancement::processing() const {
    LockGuard<Mutex> lock(_mux);
    return _num_processing;
}

SizeType WorkloadAdvancement::completed() const {
    LockGuard<Mutex> lock(_mux);
    return _num_completed;
}

SizeType WorkloadAdvancement::total() const {
    LockGuard<Mutex> lock(_mux);
    return _num_waiting + _num_processing + _num_completed;
}

Void WorkloadAdvancement::add_to_waiting(SizeType n) {
    ARIADNE_PRECONDITION(n > 0);
    LockGuard<Mutex> lock(_mux);
    _num_waiting+=n;
}

Void WorkloadAdvancement::add_to_processing(SizeType n) {
    ARIADNE_PRECONDITION(n <= _num_waiting);
    LockGuard<Mutex> lock(_mux);
    _num_waiting-=n;
    _num_processing+=n;
}

Void WorkloadAdvancement::add_to_completed(SizeType n) {
    ARIADNE_PRECONDITION(n <=_num_processing);
    LockGuard<Mutex> lock(_mux);
    _num_processing-=n;
    _num_completed+=n;
}

Double WorkloadAdvancement::completion_rate() const {
    LockGuard<Mutex> lock(_mux);
    Double total = _num_waiting + _num_processing + _num_completed;
    return ((Double)_num_completed) / total;
}

Bool WorkloadAdvancement::has_finished() const {
    LockGuard<Mutex> lock(_mux);
    return _num_processing == 0 and _num_waiting == 0;
}

} // namespace Ariadne
