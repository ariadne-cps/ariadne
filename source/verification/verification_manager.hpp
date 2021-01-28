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
#include "../concurrency/task_appraisal.hpp"
#include "../concurrency/task_appraisal_space.hpp"
#include "configuration_property_refinement_rule.hpp"

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

    template<class R> void add_safety_specification(
            TaskRunnable<R>& runnable,
            Set<TaskAppraisalConstraint<R>> const& safety_constraints,
            List<ConfigurationPropertyRefinementRule<R>> const& refinement_rules,
            Set<TimedRadiusObjective> const& objectives) const {
        auto appraisal_space = runnable.runner()->task().appraisal_space();
        auto original_constraints = appraisal_space.constraints();
        auto original_weights = appraisal_space.parameters_weights();
        TaskAppraisalSpaceBuilder<R> builder;
        for (auto constr : appraisal_space.constraints()) builder.add(constr,original_weights.get(constr.parameter()));
        for (auto constr : safety_constraints) builder.add(constr); // Use default weight 1.0
        runnable.runner()->task().set_appraisal_space(builder.build());
        runnable.runner()->task().set_configuration_refinement_rules(refinement_rules);
        runnable.runner()->task().set_configuration_refinement_objectives(objectives);
        runnable.runner()->refine_configuration_init();
    }
};

} // namespace Ariadne

#endif // ARIADNE_VERIFICATION_MANAGER_HPP
