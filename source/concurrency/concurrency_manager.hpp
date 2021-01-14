/***************************************************************************
 *            concurrency/concurrency_manager.hpp
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

/*! \file concurrency/concurrency_manager.hpp
 *  \brief Static class for managing concurrency across the library
 */

#ifndef ARIADNE_CONCURRENCY_MANAGER_HPP
#define ARIADNE_CONCURRENCY_MANAGER_HPP

#include <thread>
#include "../utility/container.hpp"
#include "../utility/pointer.hpp"

namespace Ariadne {

//! \brief Manages threads and returns runners based on concurrency availability.
class ConcurrencyManager {
  public:
    ConcurrencyManager();
    ConcurrencyManager(ConcurrencyManager const&) = delete;
    void operator=(ConcurrencyManager const&) = delete;

    static ConcurrencyManager& get_instance() {
        static ConcurrencyManager instance;
        return instance;
    }

    unsigned int maximum_concurrency() const;
    unsigned int concurrency() const;

    void set_concurrency(unsigned int value);

  private:
    const unsigned int _maximum_concurrency;
    unsigned int _concurrency;
};

} // namespace Ariadne

#endif // ARIADNE_CONCURRENCY_MANAGER_HPP
