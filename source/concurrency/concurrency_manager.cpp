/***************************************************************************
 *            concurrency/concurrency_manager.cpp
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

#include <cstdlib>
#include <time.h>

#include "../utility/macros.hpp"
#include "../concurrency/concurrency_manager.hpp"

namespace Ariadne {

ConcurrencyManager::ConcurrencyManager() : _maximum_concurrency(std::thread::hardware_concurrency()), _concurrency(1) {
    srand(time(NULL));
}

unsigned int ConcurrencyManager::maximum_concurrency() const {
    return _maximum_concurrency;
}

unsigned int ConcurrencyManager::concurrency() const {
    return _concurrency;
}

void ConcurrencyManager::set_concurrency(unsigned int value) {
    ARIADNE_PRECONDITION(value <= _maximum_concurrency and value > 0);
    std::lock_guard<std::mutex> lock(_data_mutex);
    _concurrency = value;
}

List<TaskSearchPointAppraisal> ConcurrencyManager::last_search_best_points() const {
    return _last_search_best_points;
}

void ConcurrencyManager::set_last_search_best_points(List<TaskSearchPointAppraisal> const& points) {
    std::lock_guard<std::mutex> lock(_data_mutex);
    _last_search_best_points = points;
}

}