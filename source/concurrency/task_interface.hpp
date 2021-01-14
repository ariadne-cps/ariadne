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

namespace Ariadne {

class TaskSearchPoint;
class TaskSearchPointAppraisal;
class TaskSearchSpace;
template<class I, class O> class TaskAppraisalSpace;

typedef std::chrono::microseconds DurationType;

template<class I, class O, class C>
class TaskInterface {
  public:
    typedef I InputType;
    typedef O OutputType;
    typedef C ConfigurationType;

    //! \brief The name of the task, to be used for thread naming
    virtual String name() const = 0;
    //! \brief Return the parameter space for the task
    virtual TaskSearchSpace const& search_space() const = 0;
    //! \brief Return the appraisal space for the task
    virtual TaskAppraisalSpace<I,O> const& appraisal_space() const = 0;

    //! \brief Convert a task parameter point into a configuration of values for the task, possibly using \a in for values
    virtual C to_configuration(I const& in, TaskSearchPoint const& p) const = 0;
    //! \brief The task to be performed, taking \a in as input and \a cfg as a configuration of the parameters
    virtual O run_task(I const& in, C const& cfg) const = 0;
    //! \brief Evaluate the costs of points from output and execution time, possibly using the input \a in
    virtual Set<TaskSearchPointAppraisal> appraise(Map<TaskSearchPoint,Pair<O,DurationType>> const& data, I const& in) const = 0;
};

//! \brief The base for parameter search tasks
//! \details Useful to streamline task construction
template<class I, class O, class C>
class ParameterSearchTaskBase : public TaskInterface<I,O,C> {
  public:
    typedef I InputType;
    typedef O OutputType;
  protected:
    ParameterSearchTaskBase(String const& name, TaskSearchSpace const& search_space, TaskAppraisalSpace<I,O> const& appraisal_space)
        : _name(name), _search_space(search_space), _appraisal_space(appraisal_space) {}
  public:
    String name() const override { return _name; }
    TaskSearchSpace const& search_space() const override { return _search_space; }
    TaskAppraisalSpace<I,O> const& appraisal_space() const override { return _appraisal_space; }
    Set<TaskSearchPointAppraisal> appraise(Map<TaskSearchPoint,Pair<O,DurationType>> const& data, I const& input) const override { return _appraisal_space.appraise(data,input); }
  private:
    String const _name;
    TaskSearchSpace const _search_space;
    TaskAppraisalSpace<I,O> const _appraisal_space;
};

} // namespace Ariadne

#endif // ARIADNE_TASK_INTERFACE_HPP
