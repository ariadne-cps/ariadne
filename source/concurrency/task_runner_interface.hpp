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

class TaskSearchPoint;
class TaskSearchPointCost;
class TaskSearchSpace;
template<class I, class O> class TaskIOData;

template<class I, class O, class C>
class TaskInterface {
  public:
    typedef I InputType;
    typedef O OutputType;
    typedef C ConfigurationType;

    //! \brief The name of the task, to be used for thread naming
    virtual std::string name() const = 0;
    //! \brief Return the parameter space for the task
    virtual TaskSearchSpace const& search_space() const = 0;

    //! \brief Convert a task parameter point into a configuration of values for the task, possibly using \a in for values
    virtual ConfigurationType to_configuration(InputType const& in, TaskSearchPoint const& p) const = 0;
    //! \brief The task to be performed, taking \a in as input and \a cfg as a configuration of the parameters
    virtual OutputType run_task(InputType const& in, ConfigurationType const& cfg) const = 0;
    //! \brief Evaluate the costs of points from the related data, comprising input, output and execution time
    virtual Set<TaskSearchPointCost> appraise(Map<TaskSearchPoint,TaskIOData<I,O>> const& data) const = 0;
};

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
    //! \brief Push input to the runner
    virtual void push(InputType const& input) = 0;
    //! \brief Pull output from the runner
    virtual OutputType pull() = 0;
};

//! \brief Interface for a class that supports a runnable task.
template<class T>
class TaskRunnableInterface {
  public:
    typedef T TaskType;

    //! \brief Set a new runner, useful to override the default runner
    virtual void set_runner(SharedPointer<TaskRunnerInterface<TaskType>> runner) = 0;
};

} // namespace Ariadne

#endif // ARIADNE_TASK_RUNNER_INTERFACE_HPP
