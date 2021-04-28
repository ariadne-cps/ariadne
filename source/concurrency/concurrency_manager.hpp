/***************************************************************************
 *            concurrency/concurrency_manager.hpp
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

/*! \file concurrency/concurrency_manager.hpp
 *  \brief Singleton class for managing concurrency across the library
 */

#ifndef ARIADNE_CONCURRENCY_MANAGER_HPP
#define ARIADNE_CONCURRENCY_MANAGER_HPP

#include <algorithm>
#include "utility/container.hpp"
#include "utility/pointer.hpp"
#include "concurrency/thread_pool.hpp"

namespace Ariadne {

//! \brief Manages threads based on concurrency availability.
class ConcurrencyManager {
  private:
    ConcurrencyManager();
  public:
    ConcurrencyManager(ConcurrencyManager const&) = delete;
    void operator=(ConcurrencyManager const&) = delete;

    //! \brief The singleton instance of this class
    static ConcurrencyManager& instance() {
        static ConcurrencyManager instance;
        return instance;
    }

    //! \brief Get the maximum concurrency allowed by this machine
    SizeType maximum_concurrency() const;
    //! \brief Get the preferred concurrency to be used
    //! \details A concurrency of zero is allowed, meaning that a task
    //! will be run sequentially
    SizeType concurrency() const;

    //! \brief Synchronised method for updating the preferred concurrency to be used
    void set_concurrency(SizeType value);

    //! \brief Enqueue a task for execution, returning the future handler
    //! \details The is no limits on the number of tasks to enqueue. If concurrency is zero,
    //! then the task is executed sequentially with no threads involved
    template<class F, class... AS> auto enqueue(F &&f, AS &&... args) -> Future<ResultOf<F(AS...)>>;

  private:
    const SizeType _maximum_concurrency;
    SizeType _concurrency;
    mutable Mutex _concurrency_mutex;

    ThreadPool _pool;
};

template<class F, class... AS> auto ConcurrencyManager::enqueue(F &&f, AS &&... args) -> Future<ResultOf<F(AS...)>> {
    if (_concurrency == 0) {
        using ReturnType = ResultOf<F(AS...)>;
        auto task = PackagedTask<ReturnType()>(std::bind(std::forward<F>(f), std::forward<AS>(args)...));
        Future<ReturnType> result = task.get_future();
        task();
        return result;
    } else return _pool.enqueue(f,args...);
}

} // namespace Ariadne

#endif // ARIADNE_CONCURRENCY_MANAGER_HPP
