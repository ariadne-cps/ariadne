/***************************************************************************
 *            concurrency/task_ranking_space.hpp
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

/*! \file concurrency/task_raking_space.hpp
 *  \brief Class for handling ranking of results from parameter search of a task.
 */

#ifndef ARIADNE_TASK_RANKING_SPACE_HPP
#define ARIADNE_TASK_RANKING_SPACE_HPP

#include "../algebra/vector.hpp"
#include "../concurrency/task_ranking_parameter.hpp"
#include "../concurrency/task_execution_ranking.hpp"

namespace Ariadne {

using std::min, std::max;

template<class R> class TaskRankingSpace;

template<class R>
class TaskRankingSpaceBuilder {
    typedef ScoreType WeightType;
public:
    TaskRankingSpaceBuilder() = default;

    TaskRankingSpaceBuilder& add(TaskRankingParameter<R> const& parameter, WeightType const& weight = WeightType(1.0)) {
        return add(TaskRankingConstraint<R>(parameter), weight);
    }

    TaskRankingSpaceBuilder& add(TaskRankingConstraint<R> const& constraint, WeightType const& weight = WeightType(1.0)) {
        ARIADNE_PRECONDITION(weight >= 0);
        _constraint_weights.insert(Pair<TaskRankingConstraint<R>,WeightType>(constraint, weight));
        return *this;
    }

    TaskRankingSpace<R> build() const;

  private:
    Map<TaskRankingConstraint<R>,WeightType> _constraint_weights;
};

template<class R>
class TaskRankingSpace : public WritableInterface {
    friend class TaskRankingSpaceBuilder<R>;
  public:
    typedef ScoreType WeightType;
    typedef TaskInput<R> InputType;
    typedef TaskOutput<R> OutputType;

  protected:
    TaskRankingSpace(Map<TaskRankingConstraint<R>,WeightType> const& constraint_weights)
        : _constraint_weights(constraint_weights) { }
  public:
    Map<TaskRankingConstraint<R>,WeightType> const& constraint_weights() const { return _constraint_weights; }

    TaskRankingSpace* clone() const { return new TaskRankingSpace(*this); }

    virtual OutputStream& _write(OutputStream& os) const { os << _constraint_weights; return os; }

    Bool has_critical_constraints() const {
        for (auto c : _constraint_weights.keys()) if (c.severity() == RankingConstraintSeverity::CRITICAL) return true;
        return false;
    }

    Set<TaskRankingConstraint<R>> failed_critical_constraints(InputType const& input, OutputType const& output) const {
        Set<TaskRankingConstraint<R>> result;
        for (auto c : _constraint_weights.keys()) {
            if (c.severity() == RankingConstraintSeverity::CRITICAL) {
                auto const& p = c.parameter();
                auto rank = p.rank(input,output,DurationType(0));
                auto threshold = c.threshold();
                if ((p.optimisation() == OptimisationCriterion::MINIMISE and rank > threshold) or
                    (p.optimisation() == OptimisationCriterion::MAXIMISE and rank < threshold)) {
                    result.insert(c);
                }
            }
        }
        return result;
    }

    Set<TaskExecutionRanking> rank(Map<ConfigurationSearchPoint,Pair<OutputType,DurationType>> const& data, InputType const& input) const {
        typedef TaskRankingConstraint<R> ConstraintType;
        Set<TaskExecutionRanking> result;

        // Compute the dimensions and initialise the min/max entries
        Map<ConstraintType,SizeType> dimensions;
        Map<ConstraintType,Pair<ScoreType,ScoreType>> scalar_min_max;
        Map<ConstraintType,Vector<Pair<ScoreType,ScoreType>>> vector_min_max;
        auto data_iter = data.cbegin();
        for (auto cw : _constraint_weights) {
            auto c = cw.first;
            auto p = c.parameter();
            auto dim = p.dimension(input);
            dimensions.insert(Pair<ConstraintType,SizeType>{c,dim});
            if (dim == 1) {
                auto val = p.rank(input, data_iter->second.first, data_iter->second.second);
                scalar_min_max.insert(Pair<ConstraintType,Pair<ScoreType,ScoreType>>{c, {val, val}});
            } else {
                Vector<Pair<ScoreType,ScoreType>> vals(dim);
                for (SizeType i=0; i<dim; ++i) {
                    auto val = p.rank(input, data_iter->second.first, data_iter->second.second, i);
                    vals[i] = {val,val};
                }
                vector_min_max.insert(Pair<ConstraintType,Vector<Pair<ScoreType,ScoreType>>>{c, vals});
            }
        }

        // Update the min/max on the remaining data entries
        ++data_iter;
        while (data_iter != data.cend()) {
            for (auto cw : _constraint_weights) {
                auto c = cw.first;
                auto p = c.parameter();
                auto dim = dimensions.get(c);
                if (dim == 1) {
                    auto val = p.rank(input, data_iter->second.first, data_iter->second.second);
                    scalar_min_max[c] = {min(scalar_min_max[c].first,val), max(scalar_min_max[c].second,val)};
                } else {
                    for (SizeType i=0; i<dim; ++i) {
                        auto val = p.rank(input, data_iter->second.first, data_iter->second.second, i);
                        vector_min_max[c][i] = {min(vector_min_max[c][i].first,val),max(vector_min_max[c][i].second,val)};
                    }
                }
            }
            ++data_iter;
        }

        // Compute the score
        for (auto entry : data) {
            // If max != min, compute (val-min)/(max-min), in the vector case also divide by the dimension
            // sum/subtract each weighted cost based on the minimise/maximise objective
            ScoreType score(0);
            SizeType low_errors(0), high_errors(0);
            for (auto cw : _constraint_weights) {
                auto const& c = cw.first;
                auto p = c.parameter();
                auto weight = cw.second;
                auto dim = dimensions.get(c);
                ScoreType local_score(0);
                if (dim == 1) {
                    auto max_min_diff = scalar_min_max[c].second - scalar_min_max[c].first;
                    auto rank = p.rank(input, entry.second.first, entry.second.second);
                    auto threshold = c.threshold();
                    if ((p.optimisation() == OptimisationCriterion::MINIMISE and rank > threshold) or
                        (p.optimisation() == OptimisationCriterion::MAXIMISE and rank < threshold)) {
                        if (c.severity() == RankingConstraintSeverity::PERMISSIVE) ++low_errors;
                        else if (c.severity() == RankingConstraintSeverity::CRITICAL) ++high_errors;
                    }
                    if (max_min_diff > 0) local_score = (rank-scalar_min_max[c].first)/max_min_diff;
                } else {
                    SizeType effective_dim = dim;
                    for (SizeType i=0; i<dim; ++i) {
                        auto max_min_diff = vector_min_max[c][i].second - vector_min_max[c][i].first;
                        auto rank = p.rank(input, entry.second.first, entry.second.second, i);
                        if (max_min_diff > 0) local_score = (rank - vector_min_max[c][i].first)/max_min_diff;
                        else --effective_dim;
                        if (effective_dim > 0) local_score/=effective_dim;
                    }
                }

                if (p.optimisation() == OptimisationCriterion::MAXIMISE) score += weight * local_score;
                else score -= weight*local_score;
            }
            result.insert(TaskExecutionRanking(entry.first, score, low_errors, high_errors));
        }

        return result;
    }

  private:
    Map<TaskRankingConstraint<R>,WeightType> const _constraint_weights;
};

template<class R>
TaskRankingSpace<R> TaskRankingSpaceBuilder<R>::build() const {
    ARIADNE_PRECONDITION(not _constraint_weights.empty());
    return TaskRankingSpace(_constraint_weights);
}

} // namespace Ariadne

#endif // ARIADNE_TASK_RANKING_SPACE_HPP
