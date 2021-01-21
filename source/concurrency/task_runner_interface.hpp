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
#include "../concurrency/searchable_configuration.hpp"

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

    //! \brief Activates the runner, activating the threads and consequently properly setting the verbosity level
    //! \details This is necessary since the object may need to be visible at a higher level of abstraction compared
    //! to the use. Multiple activations are possible due to multiple calls of the runner, though they will not have any effect,
    //! since a thread can be activated only once, with no exception thrown.
    virtual void activate() = 0;

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
    TaskRunnable(ConfigurationType const& configuration) : Configurable<R>(configuration) { }
    //! \brief Set a new runner, useful to override the default runner
    void set_runner(SharedPointer<TaskRunnerInterface<R>> runner) { this->_runner = runner; }
    //! \brief Get the runner
    SharedPointer<TaskRunnerInterface<R>>& runner() { return _runner; }
    SharedPointer<TaskRunnerInterface<R>> const& runner() const { return _runner; }
  private:
    SharedPointer<TaskRunnerInterface<R>> _runner;
};

} // namespace Ariadne

#endif // ARIADNE_TASK_RUNNER_INTERFACE_HPP
