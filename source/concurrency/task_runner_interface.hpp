/***************************************************************************
 *            concurrency/task_runner_interface.hpp
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

/*! \file concurrency/task_runner_interface.hpp
 *  \brief The interfaces for generic runners.
 */

#ifndef ARIADNE_TASK_RUNNER_INTERFACE_HPP
#define ARIADNE_TASK_RUNNER_INTERFACE_HPP

#include "../utility/container.hpp"
#include "../utility/pointer.hpp"
#include "../concurrency/task_interface.hpp"
#include "../concurrency/configurable.hpp"

namespace Ariadne {

//! \brief Interface for the runner of a task.
//! \details Takes the runnable class as template argument.
template<class R>
class TaskRunnerInterface {
  public:
    typedef Task<R> TaskType;
    typedef TaskInput<R> InputType;
    typedef TaskOutput<R>  OutputType;
    typedef Configuration<R> ConfigurationType;

    //! \brief Return the task
    virtual TaskType& task() = 0;
    virtual TaskType const& task() const = 0;

    //! \brief Return the configuration
    virtual ConfigurationType const& configuration() const = 0;

    //! \brief Transfer running statistics onto the ConcurrencyManager
    virtual void dump_statistics() = 0;
    //! \brief Push input
    virtual void push(InputType const& input) = 0;
    //! \brief Pull output from the runner
    virtual OutputType pull() = 0;
};

//! \brief Interface for a class that supports a runnable task.
template<class R>
class TaskRunnable : public Configurable<R> {
    friend class ConcurrencyManager;
    friend class VerificationManager;
    typedef Configuration<R> ConfigurationType;
  protected:
    TaskRunnable(ConfigurationType const& configuration);
    //! \brief Set a new runner, useful to override the default runner
    void set_runner(SharedPointer<TaskRunnerInterface<R>> runner);
    //! \brief Get the runner
    SharedPointer<TaskRunnerInterface<R>>& runner();
    SharedPointer<TaskRunnerInterface<R>> const& runner() const;
  private:
    SharedPointer<TaskRunnerInterface<R>> _runner;
};

} // namespace Ariadne

#endif // ARIADNE_TASK_RUNNER_INTERFACE_HPP
