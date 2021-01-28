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
        _ranking_constraints.insert(constraint);
        _parameters_weights.insert(Pair<TaskRankingParameter<R>,WeightType>(constraint.parameter(), weight));
        return *this;
    }

    TaskRankingSpace<R> build() const;

  private:
    Set<TaskRankingConstraint<R>> _ranking_constraints;
    Map<TaskRankingParameter<R>,WeightType> _parameters_weights;
};

template<class R>
class TaskRankingSpace : public WritableInterface {
    friend class TaskRankingSpaceBuilder<R>;
  public:
    typedef ScoreType WeightType;
    typedef TaskInput<R> InputType;
    typedef TaskOutput<R> OutputType;

  protected:
    TaskRankingSpace(Set<TaskRankingConstraint<R>> const& constraints,
                     Map<TaskRankingParameter<R>,WeightType> const& parameters_weights)
        : _constraints(constraints), _parameters_weights(parameters_weights) { }
  public:
    Set<TaskRankingConstraint<R>> const& constraints() const { return _constraints; }
    Map<TaskRankingParameter<R>,WeightType> const& parameters_weights() const { return _parameters_weights; }

    TaskRankingSpace* clone() const { return new TaskRankingSpace(*this); }

    virtual OutputStream& _write(OutputStream& os) const { os << _parameters_weights; return os; }

    Set<TaskExecutionRanking>
    rank(Map<TaskSearchPoint,Pair<OutputType,DurationType>> const& data, InputType const& input) const {
        typedef TaskRankingParameter<R> ParamType;
        Set<TaskExecutionRanking> result;

        // Compute the dimensions and initialise the min/max entries
        Map<ParamType,SizeType> dimensions;
        Map<ParamType,Pair<ScoreType,ScoreType>> scalar_min_max;
        Map<ParamType,Vector<Pair<ScoreType,ScoreType>>> vector_min_max;
        auto data_iter = data.cbegin();
        for (auto ac : _constraints) {
            auto p = ac.parameter();
            auto dim = p.dimension(input);
            dimensions.insert(Pair<ParamType,SizeType>{p,dim});
            if (dim == 1) {
                auto val = p.rank(input, data_iter->second.first, data_iter->second.second);
                scalar_min_max.insert(Pair<ParamType,Pair<ScoreType,ScoreType>>{p, {val, val}});
            } else {
                Vector<Pair<ScoreType,ScoreType>> vals(dim);
                for (SizeType i=0; i<dim; ++i) {
                    auto val = p.rank(input, data_iter->second.first, data_iter->second.second, i);
                    vals[i] = {val,val};
                }
                vector_min_max.insert(Pair<ParamType,Vector<Pair<ScoreType,ScoreType>>>{p, vals});
            }
        }

        // Update the min/max on the remaining data entries
        ++data_iter;
        while (data_iter != data.cend()) {
            for (auto ac : _constraints) {
                auto p = ac.parameter();
                auto dim = dimensions.get(p);
                if (dim == 1) {
                    auto val = p.rank(input, data_iter->second.first, data_iter->second.second);
                    scalar_min_max[p] = {min(scalar_min_max[p].first,val), max(scalar_min_max[p].second,val)};
                } else {
                    for (SizeType i=0; i<dim; ++i) {
                        auto val = p.rank(input, data_iter->second.first, data_iter->second.second, i);
                        vector_min_max[p][i] = {min(vector_min_max[p][i].first,val),max(vector_min_max[p][i].second,val)};
                    }
                }
            }
            ++data_iter;
        }

        // Compute the score
        for (auto entry : data) {
            // If max != min, compute (val-min)/(max-min), in the vector case also divide by the dimension
            // sum/subtract each weighted cost based on the minimise/maximise objective
            ScoreType c(0);
            SizeType low_errors(0), high_errors(0);
            for (auto rc : _constraints) {
                auto p = rc.parameter();
                auto weight = _parameters_weights.get(p);
                auto dim = dimensions.get(p);
                ScoreType score(0);
                if (dim == 1) {
                    auto max_min_diff = scalar_min_max[p].second - scalar_min_max[p].first;
                    auto rank = p.rank(input, entry.second.first, entry.second.second);
                    auto threshold = rc.threshold();
                    if ((p.optimisation() == RankingParameterOptimisation::MINIMISE and rank > threshold) or
                        (p.optimisation() == RankingParameterOptimisation::MAXIMISE and rank < threshold)) {
                        if (rc.severity() == RankingConstraintSeverity::PERMISSIVE) ++low_errors;
                        else if (rc.severity() == RankingConstraintSeverity::CRITICAL) ++high_errors;
                    }
                    if (max_min_diff > 0) score = (rank-scalar_min_max[p].first)/max_min_diff;
                } else {
                    SizeType effective_dim = dim;
                    for (SizeType i=0; i<dim; ++i) {
                        auto max_min_diff = vector_min_max[p][i].second - vector_min_max[p][i].first;
                        auto rank = p.rank(input, entry.second.first, entry.second.second, i);
                        if (max_min_diff > 0) score = (rank - vector_min_max[p][i].first)/max_min_diff;
                        else --effective_dim;
                        if (effective_dim > 0) c/=effective_dim;
                    }
                }

                if (p.optimisation() == RankingParameterOptimisation::MAXIMISE) c += weight * score;
                else c -= weight*score;
            }
            result.insert(TaskExecutionRanking(entry.first, c, low_errors, high_errors));
        }
        return result;
    }

  private:
    Set<TaskRankingConstraint<R>> const _constraints;
    Map<TaskRankingParameter<R>,WeightType> const _parameters_weights;
};

template<class R>
TaskRankingSpace<R> TaskRankingSpaceBuilder<R>::build() const {
    ARIADNE_PRECONDITION(not _ranking_constraints.empty());
    return TaskRankingSpace(_ranking_constraints, _parameters_weights);
}

} // namespace Ariadne

#endif // ARIADNE_TASK_RANKING_SPACE_HPP
