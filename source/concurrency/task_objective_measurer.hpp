/***************************************************************************
 *            concurrency/task_objective_measurer.hpp
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

/*! \file concurrency/task_objective_measurer.hpp
 *  \brief Classes for mea the measure of a task refinement
 */

#ifndef ARIADNE_TASK_OBJECTIVE_MEASURER_HPP
#define ARIADNE_TASK_OBJECTIVE_MEASURER_HPP

namespace Ariadne {

//! \brief A measurer for a set of task objectives related to the refinement of a property
//! \details The measurer is generic and is designed to return two measures: an error and a progress
template<class R> class TaskObjectiveMeasurerInterface {
  public:
    typedef TaskObjective<R> ObjectiveType;
    typedef TaskInput<R> InputType;
    typedef TaskOutput<R> OutputType;
    typedef double ErrorType;
    typedef double ProgressType;

    virtual Pair<ErrorType,ProgressType> get(InputType const& i, OutputType const& o, List<ObjectiveType> const& objectives) const = 0;
  protected:
    //! \brief Get the error with respect to the objectives
    //! \details The error is chosen as the maximum among the objectives, each calculated as (current-reference)/reference
    virtual ErrorType get_error(InputType const& i, OutputType const& o, List<ObjectiveType> const& objectives) const = 0;
    //! \brief Whether to discard an objective at the current input
    virtual Bool discard(InputType const& i, ObjectiveType const& obj) const = 0;
    //! \brief Get the current measure for error from input/output, where the objective may also be required to select proper data
    virtual ErrorType current_measure(InputType const& i, OutputType const& o, ObjectiveType const& obj) const = 0;
    //! \brief Get the reference measure for error from input and objective
    virtual ErrorType reference_measure(InputType const& i, ObjectiveType const& obj) const = 0;
    //! \brief Get the progress in the task with respect to the input and output
    virtual ProgressType get_progress(InputType const& i, OutputType const& o) const = 0;

    virtual ~TaskObjectiveMeasurerInterface() = default;
};

//! \brief Class to be specialised for a specific runnable, deriving from TaskObjectiveMeasurerBae
template<class R> class TaskObjectiveMeasurer;

template<class R> class TaskObjectiveMeasurerBase : public TaskObjectiveMeasurerInterface<R> {
    typedef TaskObjective<R> ObjectiveType;
    typedef TaskInput<R> InputType;
    typedef TaskOutput<R> OutputType;
    typedef typename TaskObjectiveMeasurerInterface<R>::ErrorType ErrorType;
    typedef typename TaskObjectiveMeasurerInterface<R>::ProgressType ProgressType;
  public:
    Pair<ErrorType,ProgressType> get(InputType const& i, OutputType const& o, List<ObjectiveType> const& objectives) const final {
        return make_pair(get_error(i,o,objectives),this->get_progress(i,o));
    }

    ScoreType get_error(InputType const& i, OutputType const& o, List<ObjectiveType> const& objectives) const override final {
        auto result = std::numeric_limits<double>::lowest();
        for (auto obj : objectives) {
            if (not this->discard(i,obj)) {
                auto ref = this->reference_measure(i, obj);
                auto amount = (this->current_measure(i, o, obj) - ref) / ref;
                if (amount > result) result = amount;
            }
        }
        return result;
    }
};

} // namespace Ariadne

#endif // ARIADNE_TASK_OBJECTIVE_MEASURER_HPP
