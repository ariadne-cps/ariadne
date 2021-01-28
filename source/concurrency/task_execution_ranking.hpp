/***************************************************************************
 *            concurrency/task_execution_ranking.hpp
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

/*! \file concurrency/task_execution_ranking.hpp
 *  \brief Class for defining the score for a point.
 */

#ifndef ARIADNE_TASK_EXECUTION_RANKING_HPP
#define ARIADNE_TASK_EXECUTION_RANKING_HPP

#include "../concurrency/task_search_point.hpp"
#include "../concurrency/task_ranking_parameter.hpp"

namespace Ariadne {

typedef double ScoreType;

class TaskExecutionRanking;

class CriticalRankingFailureException : public std::runtime_error {
public:
    CriticalRankingFailureException(Set<TaskExecutionRanking> const& rankings) : std::runtime_error("All ranked points have critical failures: " + to_string(rankings)) { }
};

//! \brief Enumeration for the severity of the constraint
//! \details NONE: there actually is no constraint
//!          PERMISSIVE: satisfying the constraint is only desired
//!          CRITICAL: satisfying the constraint is mandatory
enum class RankingConstraintSeverity { NONE, PERMISSIVE, CRITICAL };
inline std::ostream& operator<<(std::ostream& os, const RankingConstraintSeverity severity) {
    switch (severity) {
        case RankingConstraintSeverity::NONE: os << "NONE"; break;
        case RankingConstraintSeverity::PERMISSIVE: os << "PERMISSIVE"; break;
        case RankingConstraintSeverity::CRITICAL: os << "CRITICAL"; break;
        default: ARIADNE_FAIL_MSG("Unhandled RankingConstraintSeverity value.");
    }
    return os;
}

//! \brief Constraint for task ranking
//! \details The predicate depends on the MINIMISE/MAXIMISE character of the ranking parameter:
//! if MINIMISE, then the cost value must be lower than the threshold, if MAXIMISE it must be higher
template<class R>
class TaskRankingConstraint : public WritableInterface {
  public:
    TaskRankingConstraint(TaskRankingParameter<R> const& parameter) : _parameter(parameter), _threshold(ScoreType(0)), _severity(RankingConstraintSeverity::NONE) { }
    TaskRankingConstraint(TaskRankingParameter<R> const& parameter, ScoreType threshold, RankingConstraintSeverity severity) : _parameter(parameter), _threshold(threshold), _severity(severity) {
        ARIADNE_PRECONDITION(parameter.is_scalar() or severity == RankingConstraintSeverity::NONE); // Vector parameters do not support constraints
    }

    TaskRankingConstraint* clone() const { return new TaskRankingConstraint(*this); }
    TaskRankingParameter<R> parameter() const { return _parameter; }
    ScoreType threshold() const { return _threshold; }
    RankingConstraintSeverity severity() const { return _severity; }

    Bool operator<(TaskRankingConstraint const& c) const { return _parameter < c._parameter; }

    virtual OutputStream& _write(OutputStream& os) const {
        os << "{" << _parameter.name();
        if (_severity != RankingConstraintSeverity::NONE) {
            if (_parameter.optimisation() == OptimisationCriterion::MINIMISE) os << "<=";
            else os << ">=";
            os << _threshold << "(" << _severity << ")";
        }
        os << "}"; return os;
    }

  private:
    const TaskRankingParameter<R> _parameter;
    const ScoreType _threshold;
    const RankingConstraintSeverity _severity;
};

class TaskExecutionRanking : public WritableInterface {
public:
    TaskExecutionRanking(TaskSearchPoint const& p, ScoreType const& s, SizeType const& permissive_failures, SizeType const& critical_failures);
    TaskSearchPoint const& point() const;
    ScoreType const& score() const;
    SizeType const& permissive_failures() const;
    SizeType const& critical_failures() const;
    //! \brief Ordering is based on number of failures, followed by cost
    Bool operator<(TaskExecutionRanking const& s) const;

    virtual OutputStream& _write(OutputStream& os) const;
private:
    TaskSearchPoint _point;
    ScoreType _cost;
    SizeType _permissive_failures;
    SizeType _critical_failures;
};

} // namespace Ariadne

#endif // ARIADNE_TASK_EXECUTION_RANKING_HPP
