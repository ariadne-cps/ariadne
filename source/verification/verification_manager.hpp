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

#include "utility/container.hpp"
#include "utility/pointer.hpp"
#include "concurrency/task_runner.hpp"
#include "concurrency/task_execution_ranking.hpp"
#include "concurrency/task_ranking_space.hpp"
#include "concurrency/concurrency_manager.hpp"
#include "configuration/configuration_property_refinement.hpp"

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
                                      Set<ConfigurationPropertyRefinement<R>> const& refinement_targets) const
    {
        auto appraisal_space = runnable.runner()->task().ranking_space();
        auto original_constraints = appraisal_space.constraints();
        auto original_weights = appraisal_space.parameters_weights();
        TaskRankingSpaceBuilder<R> builder;
        for (auto constr : appraisal_space.constraints()) builder.add(constr,original_weights.get(constr.parameter()));
        for (auto constr : safety_constraints) builder.add(constr,1.0/safety_constraints.size()); // Use a total weight of 1.5 distributed along safety constraints
        auto search_space_dimension = runnable.searchable_configuration().search_space().dimension();
        auto num_refinements = refinement_targets.size();
        ARIADNE_ASSERT_MSG(num_refinements <= search_space_dimension, "The number of properties to refine is greater than the number of non-single properties.");
        Configuration<R> cfg = runnable.configuration();
        for (auto target : refinement_targets)
            cfg.properties().get(target.path().first())->refine_init(target.path().subpath());

        ConcurrencyManager::instance().choose_runner_for(runnable, cfg); // Re-chooses the runner
        runnable.runner()->task().set_ranking_space(builder.build());
        runnable.runner()->task().set_configuration_refinements(refinement_targets);
    }
};

} // namespace Ariadne

#endif // ARIADNE_VERIFICATION_MANAGER_HPP
