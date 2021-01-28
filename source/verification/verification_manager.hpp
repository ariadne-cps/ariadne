/***************************************************************************
 *            verification/verification_manager.hpp
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

/*! \file verification/verification_manager.hpp
 *  \brief Singleton class for managing verification processes
 */

#ifndef ARIADNE_VERIFICATION_MANAGER_HPP
#define ARIADNE_VERIFICATION_MANAGER_HPP

#include "../utility/container.hpp"
#include "../utility/pointer.hpp"
#include "../concurrency/task_runner.hpp"
#include "../concurrency/task_execution_ranking.hpp"
#include "../concurrency/task_ranking_space.hpp"
#include "configuration/configuration_property_refinement_target.hpp"

namespace Ariadne {

class VerificationManager {
  private:
    VerificationManager();
  public:
    VerificationManager(VerificationManager const&) = delete;
    void operator=(VerificationManager const&) = delete;

    static VerificationManager& instance() {
        static VerificationManager instance;
        return instance;
    }

    template<class R> void add_safety_specification(TaskRunnable<R>& runnable,
                                      Set<TaskRankingConstraint<R>> const& safety_constraints,
                                      List<ConfigurationPropertyRefinementTarget<R>> const& refinement_targets) const
    {
        auto appraisal_space = runnable.runner()->task().ranking_space();
        auto original_constraints = appraisal_space.constraints();
        auto original_weights = appraisal_space.parameters_weights();
        TaskRankingSpaceBuilder<R> builder;
        for (auto constr : appraisal_space.constraints()) builder.add(constr,original_weights.get(constr.parameter()));
        for (auto constr : safety_constraints) builder.add(constr); // Use default weight 1.0
        runnable.runner()->task().set_ranking_space(builder.build());
        runnable.runner()->task().set_configuration_refinement_targets(refinement_targets);
        runnable.runner()->refine_configuration_init();
    }
};

} // namespace Ariadne

#endif // ARIADNE_VERIFICATION_MANAGER_HPP
