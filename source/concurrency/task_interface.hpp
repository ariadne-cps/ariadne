/***************************************************************************
 *            concurrency/task_interface.hpp
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

/*! \file concurrency/task_interface.hpp
 *  \brief The interface for tasks.
 */

#ifndef ARIADNE_TASK_INTERFACE_HPP
#define ARIADNE_TASK_INTERFACE_HPP

#include "../utility/container.hpp"
#include "../utility/pointer.hpp"
#include "../utility/string.hpp"

namespace Ariadne {

class TaskSearchPoint;
class TaskSearchPointAppraisal;
class TaskSearchSpace;
template<class R> class TaskAppraisalSpace;

typedef std::chrono::microseconds DurationType;

template<class R> struct TaskInput;
template<class R> struct TaskOutput;
template<class R> class Task;
template<class R> class Configuration;
template<class R> class ConfigurationRefinementRule;

template<class R>
class TaskInterface {
  public:
    typedef TaskInput<R> InputType;
    typedef TaskOutput<R> OutputType;
    typedef Configuration<R> ConfigurationType;

    //! \brief The name of the task, to be used for thread naming
    virtual String name() const = 0;
    //! \brief Return the appraisal space for the task
    virtual TaskAppraisalSpace<R> const& appraisal_space() const = 0;
    //! \brief Set the appraisal space for the task
    virtual Void set_appraisal_space(TaskAppraisalSpace<R> const& space) = 0;
    //! \brief Return the configuration refinement rules
    virtual Set<ConfigurationRefinementRule<R>> const& configuration_refinement_rules() const = 0;
    //! \brief Set the configuration refinement rules for the task
    virtual Void set_configuration_refinement_rules(Set<ConfigurationRefinementRule<R>> const& rules) = 0;

    //! \brief The task to be performed, taking \a in as input and \a cfg as a configuration of the parameters
    virtual OutputType run_task(InputType const& in, ConfigurationType const& cfg) const = 0;
    //! \brief Evaluate the costs of points from output and execution time, possibly using the input \a in
    virtual Set<TaskSearchPointAppraisal> appraise(Map<TaskSearchPoint,Pair<OutputType,DurationType>> const& data, InputType const& in) const = 0;
};

} // namespace Ariadne

#endif // ARIADNE_TASK_INTERFACE_HPP
