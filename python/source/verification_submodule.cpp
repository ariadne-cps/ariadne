/***************************************************************************
 *            verification_submodule.cpp
 *
 *  Copyright  2008-24  Luca Geretti
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

#include "pybind11.hpp"
#include "utilities.hpp"

#include "geometry/paving_interface.hpp"
#include "geometry/grid_paving.hpp"
#include "io/figure.hpp"

#include "verification/reach_avoid.hpp"

using namespace Ariadne;

Void export_reach_avoid(pybind11::module& module)
{
    auto const& reference_internal = pybind11::return_value_policy::reference_internal;

    pybind11::class_<ReachAvoid> reach_avoid_class(module,"ReachAvoid");
    reach_avoid_class.def(pybind11::init<String,EffectiveVectorMultivariateFunction const&,Grid const&,BoundsBoxType const&,Grid const&,BoundsBoxType const&,SizeType,ExactDouble,ProbabilityType>());
    reach_avoid_class.def("add_obstacle", &ReachAvoid::add_obstacle,reference_internal);
    reach_avoid_class.def("add_goal", &ReachAvoid::add_goal,reference_internal);

    reach_avoid_class.def("state_size", &ReachAvoid::state_size);
    reach_avoid_class.def("control_size", &ReachAvoid::control_size);
    reach_avoid_class.def("obstacles_size", &ReachAvoid::obstacles_size);
    reach_avoid_class.def("goals_size", &ReachAvoid::goals_size);
    reach_avoid_class.def("unverified_size", &ReachAvoid::unverified_size);

    reach_avoid_class.def("num_sources", &ReachAvoid::num_sources);
    reach_avoid_class.def("num_destinations", &ReachAvoid::num_destinations);
    reach_avoid_class.def("num_transitions", &ReachAvoid::num_transitions);

    reach_avoid_class.def("plot", (void(ReachAvoid::*)(SizeType,SizeType)const) &ReachAvoid::plot);
    reach_avoid_class.def("plot", (void(ReachAvoid::*)()const) &ReachAvoid::plot);

    reach_avoid_class.def("print_goals", &ReachAvoid::print_goals);
    reach_avoid_class.def("print_obstacles", &ReachAvoid::print_obstacles);
    reach_avoid_class.def("print_graph", &ReachAvoid::print_graph);

    reach_avoid_class.def("compute_reachability_graph", &ReachAvoid::compute_reachability_graph);
    reach_avoid_class.def("refine_to_safety_graph", &ReachAvoid::refine_to_safety_graph);
    reach_avoid_class.def("refine_to_goal_reachable_graph", &ReachAvoid::refine_to_goal_reachable_graph);
    reach_avoid_class.def("update_unverified", &ReachAvoid::update_unverified);

    reach_avoid_class.def("unverified_percentage", &ReachAvoid::unverified_percentage);
}

Void verification_submodule(pybind11::module& module)
{
    export_reach_avoid(module);
}

