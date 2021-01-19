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

#include "../algebra/vector.hpp"
#include "../concurrency/task_appraisal_parameter.hpp"
#include "../concurrency/task_appraisal.hpp"

namespace Ariadne {

using std::min, std::max;

template<class R> class TaskAppraisalSpace;

template<class R>
class TaskAppraisalSpaceBuilder {
    typedef CostType WeightType;
public:
    TaskAppraisalSpaceBuilder() = default;

    TaskAppraisalSpaceBuilder& add(TaskAppraisalParameter<R> const& parameter, WeightType const& weight = WeightType(1.0)) {
        return add(TaskAppraisalConstraint<R>(parameter),weight);
    }

    TaskAppraisalSpaceBuilder& add(TaskAppraisalConstraint<R> const& constraint, WeightType const& weight = WeightType(1.0)) {
        ARIADNE_PRECONDITION(weight >= 0);
        _appraisal_constraints.insert(constraint);
        _parameters_weights.insert(Pair<TaskAppraisalParameter<R>,WeightType>(constraint.parameter(),weight));
        return *this;
    }

    TaskAppraisalSpace<R> build() const;

  private:
    Set<TaskAppraisalConstraint<R>> _appraisal_constraints;
    Map<TaskAppraisalParameter<R>,WeightType> _parameters_weights;
};

template<class R>
class TaskAppraisalSpace : public WritableInterface {
    friend class TaskAppraisalSpaceBuilder<R>;
  public:
    typedef CostType WeightType;
    typedef TaskInput<R> InputType;
    typedef TaskOutput<R> OutputType;

  protected:
    TaskAppraisalSpace(Set<TaskAppraisalConstraint<R>> const& appraisal_constraints,
                       Map<TaskAppraisalParameter<R>,WeightType> const& parameters_weights)
        : _appraisal_constraints(appraisal_constraints), _parameters_weights(parameters_weights) { }
  public:
    Set<TaskAppraisalConstraint<R>> const& constraints() const { return _appraisal_constraints; }
    Map<TaskAppraisalParameter<R>,WeightType> const& parameters_weights() const { return _parameters_weights; }

    TaskAppraisalSpace* clone() const { return new TaskAppraisalSpace(*this); }

    virtual OutputStream& _write(OutputStream& os) const { os << _parameters_weights; return os; }

    Set<TaskSearchPointAppraisal>
    appraise(Map<TaskSearchPoint,Pair<OutputType,DurationType>> const& data, InputType const& input) const {
        typedef TaskAppraisalParameter<R> ParamType;
        Set<TaskSearchPointAppraisal> result;

        // Compute the dimensions and initialise the min/max entries
        Map<ParamType,SizeType> dimensions;
        Map<ParamType,Pair<CostType,CostType>> scalar_min_max;
        Map<ParamType,Vector<Pair<CostType,CostType>>> vector_min_max;
        auto data_iter = data.cbegin();
        for (auto ac : _appraisal_constraints) {
            auto p = ac.parameter();
            auto dim = p.dimension(input);
            dimensions.insert(Pair<ParamType,SizeType>{p,dim});
            if (dim == 1) {
                auto val = p.appraise(input,data_iter->second.first,data_iter->second.second);
                scalar_min_max.insert(Pair<ParamType,Pair<CostType,CostType>>{p,{val,val}});
            } else {
                Vector<Pair<CostType,CostType>> vals(dim);
                for (SizeType i=0; i<dim; ++i) {
                    auto val = p.appraise(input,data_iter->second.first,data_iter->second.second,i);
                    vals[i] = {val,val};
                }
                vector_min_max.insert(Pair<ParamType,Vector<Pair<CostType,CostType>>>{p,vals});
            }
        }

        // Update the min/max on the remaining data entries
        ++data_iter;
        while (data_iter != data.cend()) {
            for (auto ac : _appraisal_constraints) {
                auto p = ac.parameter();
                auto dim = dimensions.get(p);
                if (dim == 1) {
                    auto val = p.appraise(input,data_iter->second.first,data_iter->second.second);
                    scalar_min_max[p] = {min(scalar_min_max[p].first,val), max(scalar_min_max[p].second,val)};
                } else {
                    for (SizeType i=0; i<dim; ++i) {
                        auto val = p.appraise(input,data_iter->second.first,data_iter->second.second,i);
                        vector_min_max[p][i] = {min(vector_min_max[p][i].first,val),max(vector_min_max[p][i].second,val)};
                    }
                }
            }
            ++data_iter;
        }

        // Compute the cost
        for (auto entry : data) {
            // If max != min, compute (val-min)/(max-min), in the vector case also divide by the dimension
            // sum/subtract each weighted cost based on the minimise/maximise objective
            CostType c(0);
            SizeType low_errors(0), high_errors(0);
            for (auto ac : _appraisal_constraints) {
                auto p = ac.parameter();
                auto weight = _parameters_weights.get(p);
                auto dim = dimensions.get(p);
                CostType val(0);
                if (dim == 1) {
                    auto max_min_diff = scalar_min_max[p].second - scalar_min_max[p].first;
                    auto appr = p.appraise(input,entry.second.first,entry.second.second);
                    auto threshold = ac.threshold();
                    if ((p.optimisation() == TaskAppraisalParameterOptimisation::MINIMISE and appr > threshold) or
                        (p.optimisation() == TaskAppraisalParameterOptimisation::MAXIMISE and appr < threshold)) {
                        if (ac.severity() == AppraisalConstraintSeverity::PERMISSIVE) ++low_errors;
                        else if (ac.severity() == AppraisalConstraintSeverity::CRITICAL) ++high_errors;
                    }
                    if (max_min_diff > 0) val = (appr - scalar_min_max[p].first)/max_min_diff;
                } else {
                    SizeType effective_dim = dim;
                    for (SizeType i=0; i<dim; ++i) {
                        auto max_min_diff = vector_min_max[p][i].second - vector_min_max[p][i].first;
                        auto appr = p.appraise(input,entry.second.first,entry.second.second,i);
                        if (max_min_diff > 0) val = (appr - vector_min_max[p][i].first)/max_min_diff;
                        else --effective_dim;
                        if (effective_dim > 0) c/=effective_dim;
                    }
                }

                if (p.optimisation() == TaskAppraisalParameterOptimisation::MINIMISE) c += weight*val;
                else c -= weight*val;
            }
            result.insert(TaskSearchPointAppraisal(entry.first,c,low_errors,high_errors));
        }
        return result;
    }

  private:
    Set<TaskAppraisalConstraint<R>> const _appraisal_constraints;
    Map<TaskAppraisalParameter<R>,WeightType> const _parameters_weights;
};

template<class R>
TaskAppraisalSpace<R> TaskAppraisalSpaceBuilder<R>::build() const {
    ARIADNE_PRECONDITION(not _appraisal_constraints.empty());
    return TaskAppraisalSpace(_appraisal_constraints,_parameters_weights);
}

} // namespace Ariadne

#endif // ARIADNE_TASK_APPRAISAL_SPACE_HPP
