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

#include <thread>
#include "../utility/container.hpp"
#include "../utility/pointer.hpp"
#include "../concurrency/task_runner.hpp"
#include "../concurrency/task_appraisal.hpp"
#include "../concurrency/task_appraisal_space.hpp"

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

    template<class T> void verify_runnable(TaskRunnable<T>& runnable, Set<TaskAppraisalConstraint<typename T::InputType,typename T::OutputType>> const& specification) const {
        typedef typename T::InputType I;
        typedef typename T::OutputType O;
        auto appraisal_space = runnable.runner()->task().appraisal_space();
        auto original_constraints = appraisal_space.constraints();
        auto original_weights = appraisal_space.parameters_weights();
        TaskAppraisalSpaceBuilder<I,O> builder;
        for (auto constr : appraisal_space.constraints()) {
            builder.add(constr,original_weights.get(constr.parameter()));
        }
        for (auto constr : specification) {
            builder.add(constr,0.0);
        }
        runnable.runner()->task().set_appraisal_space(builder.build());
        ARIADNE_LOG_PRINTLN("new_appraisal_space = " << runnable.runner()->task().appraisal_space());
    }
};

} // namespace Ariadne

#endif // ARIADNE_VERIFICATION_MANAGER_HPP
