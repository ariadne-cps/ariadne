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
    //static const size_t NUM_RUNS_PER_SYSTEM = 10;

    for (auto const& s : specs) {
        CONCLOG_PRINTLN(s.name)

        CONCLOG_PRINTLN_AT(1,"Generating the constraints...")
        auto constraints_prescriptions = generate_ellipsoidal_hyperbolic_constraints(NUM_CONSTRAINTS,s,configuration);
        CONCLOG_PRINTLN_AT(1,"Expected frequencies: " << frequencies(constraints_prescriptions))

        CONCLOG_PRINTLN_AT(1,"Running...")
        auto result = constrained_evolution(s.dynamics,s.initial_set,s.evolution_time,configuration,constraints(constraints_prescriptions));

        CONCLOG_PRINTLN_VAR(result.satisfaction().prescription_ratios())
        CONCLOG_PRINTLN_VAR(result.satisfaction().success_ratios())
        CONCLOG_PRINTLN_VAR(result.satisfaction().global_success_ratio())
    }

}
