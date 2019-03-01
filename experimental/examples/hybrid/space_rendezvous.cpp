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

#include "space_rendezvous.hpp"
#include <cstring>

namespace Ariadne {

using std::cout;

void verify_space_rendezvous() {
    auto problem = make_space_rendezvous_problem();

    auto initial_set = problem.initial_set;
    auto system = problem.system;
    HybridConstraintSet safe_set = problem.safe_set;

    cout << "initial_set=" << initial_set << "\n";
    cout << "system=" << system << "\n";
    cout << "safe_set=" << safe_set << "\n";

    cout << "state_space=" << system.state_space() << "\n";
    DiscreteLocation initial_location = initial_set.location();
    cout << "initial_location=" << initial_location << "\n";
    RealSpace initial_space = system.state_space()[initial_location];
    cout << "initial_space=" << initial_space << "\n";
    BoundedConstraintSet initial_constraint_set = initial_set.euclidean_set(initial_location,initial_space);
    cout << "initial_constraint_set=" << initial_constraint_set << "\n";

    MaximumError max_err=0.01;
    TaylorSeriesIntegrator integrator(max_err,Order(3u));

    GeneralHybridEvolver evolver(system);
    evolver.set_integrator(integrator);
    evolver.configuration().set_maximum_step_size(8);
    evolver.verbosity=1;

    cout << "\nComputing orbit...\n";
    HybridTime evolution_time(200.0,3);
    auto orbit=evolver.orbit(initial_set,evolution_time,Semantics::UPPER);

    RealVariable t("t"), x("x"), y("y");

    cout << "\nDrawing orbit...\n";
    plot("space_rendezvous_x_y",{-1000<=x<=200,-450<=y<=0},Colour(.5,.0,.5),orbit.reach());
}


void test() {

    verify_space_rendezvous();
}

} // namespace Ariadne

int main() {
    Ariadne::test();
}
