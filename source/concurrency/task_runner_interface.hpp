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

namespace Ariadne {

//! \brief Interface for the runner of a task.
//! \details This will usually be first implemented into an abstract base class.
template<class T>
class TaskRunnerInterface {
  public:
    typedef typename T::InputType InputType;
    typedef typename T::OutputType OutputType;
    typedef typename T::ConfigurationType ConfigurationType;

    //! \brief Activates the runner, activating the threads and consequently properly setting the verbosity level
    //! \details This is necessary since the object may need to be visible at a higher level of abstraction compared
    //! to the use. Multiple activations are possible due to multiple calls of the runner, though they will not have any effect,
    //! since a thread can be activated only once, with no exception thrown.
    virtual void activate() = 0;

    //! \brief Return the task
    virtual T& task() = 0;
    virtual T const& task() const = 0;

    //! \brief Transfer running statistics onto the ConcurrencyManager
    virtual void dump_statistics() = 0;
    //! \brief Push input to the runner
    virtual void push(InputType const& input) = 0;
    //! \brief Pull output from the runner
    virtual OutputType pull() = 0;
};

//! \brief Interface for a class that supports a runnable task.
template<class T>
class TaskRunnable {
    friend class ConcurrencyManager;
    friend class VerificationManager;
  public:
    typedef T TaskType;
  protected:
    //! \brief Set a new runner, useful to override the default runner
    void set_runner(SharedPointer<TaskRunnerInterface<TaskType>> runner) { this->_runner = runner; }
    //! \brief Get the runner
    SharedPointer<TaskRunnerInterface<TaskType>>& runner() { return _runner; }
    SharedPointer<TaskRunnerInterface<TaskType>> const& runner() const { return _runner; }
  private:
    SharedPointer<TaskRunnerInterface<TaskType>> _runner;
};

} // namespace Ariadne

#endif // ARIADNE_TASK_RUNNER_INTERFACE_HPP
