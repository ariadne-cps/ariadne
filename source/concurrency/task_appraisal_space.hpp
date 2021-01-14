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

using std::min, std::max;

template<class I, class O>
class TaskAppraisalSpace : public WritableInterface {
  public:
    typedef I InputType;
    typedef O OutputType;

    TaskAppraisalSpace(Set<TaskAppraisalParameter<I,O>> const& parameters) : _parameters(parameters) {
        ARIADNE_PRECONDITION(not _parameters.empty());
    }

    Set<TaskAppraisalParameter<I,O>> const& parameters() const { return _parameters; }

    TaskAppraisalSpace* clone() const { return new TaskAppraisalSpace(*this); }

    virtual OutputStream& _write(OutputStream& os) const { os << _parameters; return os; }

    Set<TaskSearchPointAppraisal>
    appraise(Map<TaskSearchPoint,Pair<O,DurationType>> const& data, I const& input) const {
        typedef TaskAppraisalParameter<I,O> ParamType;
        Set<TaskSearchPointAppraisal> result;

        // Compute the dimensions and initialise the min/max entries
        Map<ParamType,SizeType> dimensions;
        Map<ParamType,Pair<CostType,CostType>> scalar_min_max;
        Map<ParamType,Vector<Pair<CostType,CostType>>> vector_min_max;
        auto data_iter = data.cbegin();
        for (auto p : _parameters) {
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
            for (auto p : _parameters) {
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
            for (auto p : _parameters) {
                auto dim = dimensions.get(p);
                CostType val(0);
                if (dim == 1) {
                    auto max_min_diff = scalar_min_max[p].second - scalar_min_max[p].first;
                    if (max_min_diff > 0) val = (p.appraise(input,entry.second.first,entry.second.second) - scalar_min_max[p].first)/max_min_diff;
                } else {
                    SizeType effective_dim = dim;
                    for (SizeType i=0; i<dim; ++i) {
                        auto max_min_diff = vector_min_max[p][i].second - vector_min_max[p][i].first;
                        if (max_min_diff > 0) val = (p.appraise(input,entry.second.first,entry.second.second,i) - vector_min_max[p][i].first)/max_min_diff;
                        else --effective_dim;
                        if (effective_dim > 0) c/=effective_dim;
                    }
                }
                if (p.optimisation() == TaskAppraisalParameterOptimisation::MINIMISE) c += p.weight()*val;
                else c -= p.weight()*val;
            }
            result.insert(TaskSearchPointAppraisal(entry.first,c,0)); // (temporary) set all threshold failures to zeros
        }
        return result;
    }

  private:
    const Set<TaskAppraisalParameter<I,O>> _parameters;
};

} // namespace Ariadne

#endif // ARIADNE_TASK_APPRAISAL_SPACE_HPP
