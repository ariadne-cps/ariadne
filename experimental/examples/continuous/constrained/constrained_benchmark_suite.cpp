/***************************************************************************
 *            constrained_benchmark_suite.cpp
 *
 *  Copyright  2023  Luca Geretti
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

#include "ariadne_main.hpp"
#include "jet-engine_c.hpp"
#include "vanderpol_c.hpp"
#include "laub-loomis_c.hpp"
#include "lotka-volterra_c.hpp"
#include "higgins-selkov_c.hpp"
#include "chemical-reactor_c.hpp"
#include "pi-controller_c.hpp"
#include "rossler-attractor_c.hpp"

void ariadne_main() {
    List<SystemSpecification> specs;
    specs.push_back(JET_c());
    //specs.push_back(VDP_c());
    //specs.push_back(LAU_c());
    //specs.push_back(LOT_c());
    //specs.push_back(HIG_c());
    //specs.push_back(CHE_c());
    //specs.push_back(PIC_c());
    //specs.push_back(ROS_c());

    auto configuration = get_configuration();

    static const size_t NUM_CONSTRAINTS = 100;
    static const size_t NUM_EXECUTIONS_PER_SYSTEM = 1;

    for (auto const& s : specs) {
        CONCLOG_PRINTLN(s.name)
        for (size_t i=0; i<NUM_EXECUTIONS_PER_SYSTEM; ++i) {
            CONCLOG_PRINTLN_AT(1,"Execution #" << i)
            CONCLOG_PRINTLN_AT(2,"Generating the constraints...")
            auto constraints_prescriptions = generate_ellipsoidal_hyperbolic_constraints(NUM_CONSTRAINTS,s,configuration);
            CONCLOG_PRINTLN_AT(2,"Expected frequencies: " << frequencies(constraints_prescriptions))

            CONCLOG_PRINTLN_AT(2,"Constrained running...")
            CONCLOG_RUN_AT(1,auto c_result = constrained_evolution(s.dynamics,s.initial_set,s.evolution_time,configuration,constraints(constraints_prescriptions)))

            for (auto const& ss : c_result.satisfaction().snapshots()) {
                CONCLOG_PRINTLN_AT(2,ss.time() << ": " << ss.success_ratios() << " (" << round(ss.global_success_ratio()*100) << "%)")
            }
            CONCLOG_PRINTLN_AT(1,"Constrained cost: " << c_result.satisfaction().cost())

            CONCLOG_PRINTLN_AT(1,"Terminated in " << c_result.satisfaction().execution_time() << " s")

            CONCLOG_PRINTLN_AT(2,"Unconstrained running...")
            CONCLOG_RUN_AT(1,auto u_result = unconstrained_evolution(s.dynamics,s.initial_set,s.evolution_time,configuration,constraints(constraints_prescriptions),c_result.satisfaction().execution_time()))

            for (auto const& ss : u_result.satisfaction().snapshots()) {
                CONCLOG_PRINTLN_AT(2,ss.time() << ": " << ss.success_ratios() << " (" << round(ss.global_success_ratio()*100) << "%)")
            }
            CONCLOG_PRINTLN_AT(1,"Unconstrained cost: " << u_result.satisfaction().cost())
        }
    }

}
