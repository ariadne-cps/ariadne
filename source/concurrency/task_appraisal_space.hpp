/***************************************************************************
 *            concurrency/task_appraisal_space.hpp
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

/*! \file concurrency/task_appraisal_space.hpp
 *  \brief Class for handling appraisal of results from parameter search of a task.
 */

#ifndef ARIADNE_TASK_APPRAISAL_SPACE_HPP
#define ARIADNE_TASK_APPRAISAL_SPACE_HPP

#include "concurrency/task_appraisal_parameter.hpp"

namespace Ariadne {

template<class I, class O>
class TaskAppraisalSpace : public WritableInterface {
  public:
    typedef I InputType;
    typedef O OutputType;

    TaskAppraisalSpace(Set<TaskAppraisalParameter<InputType,OutputType>> const& parameters) : _parameters(parameters) {
        ARIADNE_PRECONDITION(not _parameters.empty());
    }

    Set<TaskAppraisalParameter<InputType,OutputType>> const& parameters() const { return _parameters; }

    TaskAppraisalSpace* clone() const { return new TaskAppraisalSpace(*this); }

    virtual OutputStream& _write(OutputStream& os) const { os << _parameters; return os; }

  private:
    const Set<TaskAppraisalParameter<InputType,OutputType>> _parameters;
};

} // namespace Ariadne

#endif // ARIADNE_TASK_APPRAISAL_SPACE_HPP
