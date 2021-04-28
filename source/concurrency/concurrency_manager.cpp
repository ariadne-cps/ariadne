/***************************************************************************
 *            concurrency/concurrency_manager.cpp
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

#include "utility/macros.hpp"
#include "concurrency/concurrency_manager.hpp"

namespace Ariadne {

ConcurrencyManager::ConcurrencyManager() : _maximum_concurrency(std::thread::hardware_concurrency()), _concurrency(0) {}

SizeType ConcurrencyManager::maximum_concurrency() const {
    return _maximum_concurrency;
}

SizeType ConcurrencyManager::concurrency() const {
    return _concurrency;
}

void ConcurrencyManager::set_concurrency(SizeType value) {
    ARIADNE_PRECONDITION(value <= _maximum_concurrency);
    LockGuard<Mutex> lock(_concurrency_mutex);
    _concurrency = value;
}

}