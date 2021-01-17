/***************************************************************************
 *            concurrency/task_appraisal.hpp
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

/*! \file concurrency/task_appraisal.hpp
 *  \brief Class for defining task appraisal classes.
 */

#ifndef ARIADNE_TASK_APPRAISAL_HPP
#define ARIADNE_TASK_APPRAISAL_HPP

#include "../concurrency/task_search_point.hpp"
#include "../concurrency/task_appraisal_parameter.hpp"

namespace Ariadne {

typedef double CostType;

//! \brief Enumeration for the severity of the constraint
//! \details NONE: there actually is no constraint
//!          LOW: satisfying the constraint is only desired
//!          HIGH: satisfying the constraint is mandatory
enum class AppraisalConstraintSeverity { NONE, LOW, HIGH };
inline std::ostream& operator<<(std::ostream& os, const AppraisalConstraintSeverity severity) {
    switch (severity) {
        case AppraisalConstraintSeverity::NONE: os << "NONE"; break;
        case AppraisalConstraintSeverity::LOW: os << "LOW"; break;
        case AppraisalConstraintSeverity::HIGH: os << "HIGH"; break;
        default: ARIADNE_FAIL_MSG("Unhandled AppraisalConstraintSeverity value.");
    }
    return os;
}

//! \brief Threshold for task appraisal
//! \details The predicate depends on the MINIMISE/MAXIMISE character of the appraisal parameter:
//! if MINIMISE, then the appraisal value must be lower than the threshold, if MAXIMISE it must be higher
template<class I, class O>
class TaskAppraisalConstraint : public WritableInterface {
  public:
    TaskAppraisalConstraint(TaskAppraisalParameter<I,O> const& parameter) : _parameter(parameter), _threshold(CostType(0)), _severity(AppraisalConstraintSeverity::NONE) { }
    TaskAppraisalConstraint(TaskAppraisalParameter<I,O> const& parameter, CostType threshold, AppraisalConstraintSeverity severity) : _parameter(parameter), _threshold(threshold), _severity(severity) {
        ARIADNE_PRECONDITION(parameter.is_scalar() or severity == AppraisalConstraintSeverity::NONE); // Vector parameters do not support constraints
    }

    TaskAppraisalConstraint* clone() const { return new TaskAppraisalConstraint(*this); }
    TaskAppraisalParameter<I,O> parameter() const { return _parameter; }
    CostType threshold() const { return _threshold; }
    AppraisalConstraintSeverity severity() const { return _severity; }

    Bool operator<(TaskAppraisalConstraint const& c) const { return _parameter < c._parameter; }

    virtual OutputStream& _write(OutputStream& os) const {
        os << "{" << _parameter.name();
        if (_severity != AppraisalConstraintSeverity::NONE) {
            if (_parameter.optimisation() == TaskAppraisalParameterOptimisation::MINIMISE) os << "<=";
            else os << ">=";
            os << _threshold << "(" << _severity << ")";
        }
        os << "}"; return os;
    }

  private:
    const TaskAppraisalParameter<I,O> _parameter;
    const CostType _threshold;
    const AppraisalConstraintSeverity _severity;
};

class TaskSearchPointAppraisal : public WritableInterface {
public:
    TaskSearchPointAppraisal(TaskSearchPoint const& p, CostType const& s, SizeType const& low_failures, SizeType const& high_failures);
    TaskSearchPoint const& point() const;
    CostType const& cost() const;
    SizeType const& low_failures() const;
    SizeType const& high_failures() const;
    //! \brief Ordering is based on number of failures, followed by cost
    Bool operator<(TaskSearchPointAppraisal const& s) const;

    virtual OutputStream& _write(OutputStream& os) const;
private:
    TaskSearchPoint _point;
    CostType _cost;
    SizeType _low_failures;
    SizeType _high_failures;
};

} // namespace Ariadne

#endif // ARIADNE_TASK_APPRAISAL_HPP
