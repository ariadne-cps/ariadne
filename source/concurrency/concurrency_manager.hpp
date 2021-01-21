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
 *  \brief Singleton class for managing concurrency across the library
 */

#ifndef ARIADNE_CONCURRENCY_MANAGER_HPP
#define ARIADNE_CONCURRENCY_MANAGER_HPP

#include <thread>
#include <algorithm>
#include "../utility/container.hpp"
#include "../utility/pointer.hpp"
#include "../concurrency/task_runner.hpp"
#include "../concurrency/task_appraisal.hpp"

namespace Ariadne {

//! \brief Manages threads and sets runners based on concurrency availability.
class ConcurrencyManager {
  private:
    ConcurrencyManager();
  public:
    ConcurrencyManager(ConcurrencyManager const&) = delete;
    void operator=(ConcurrencyManager const&) = delete;

    static ConcurrencyManager& instance() {
        static ConcurrencyManager instance;
        return instance;
    }

    template<class T> void set_runner(TaskRunnable<T>& runnable) const {
        auto const& cfg = runnable.configuration();
        SharedPointer<TaskRunnerInterface<T>> runner;
        if (_concurrency > 1 and not cfg.is_singleton())
            runner.reset(new ParameterSearchRunner<T>(cfg,std::min(_concurrency,cfg.search_space().total_points())));
        else
            runner.reset(new SequentialRunner<T>(cfg));
        runnable.set_runner(runner);
    }

    SizeType maximum_concurrency() const;
    SizeType concurrency() const;

    void set_concurrency(SizeType value);

    List<TaskSearchPointAppraisal> last_search_best_points() const;
    void set_last_search_best_points(List<TaskSearchPointAppraisal> const& points);

  private:
    const SizeType _maximum_concurrency;
    SizeType _concurrency;
    std::mutex _data_mutex;
    List<TaskSearchPointAppraisal> _last_search_best_points;
};

} // namespace Ariadne

#endif // ARIADNE_CONCURRENCY_MANAGER_HPP
