/***************************************************************************
 *            test_hybrid_reachability.cpp
 *
 *  Copyright  2006-20  Luca Geretti
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

#include <fstream>
#include <iostream>

#include "config.hpp"
#include "utility/tuple.hpp"
#include "numeric/rational.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "function/function.hpp"
#include "geometry/box.hpp"
#include "set/list_set.hpp"
#include "integration/integrator.hpp"
#include "dynamics/orbit.hpp"
#include "system/hybrid_automata.hpp"
#include "evolution/hybrid_time.hpp"
#include "hybrid/hybrid_paving.hpp"
#include "evolution/hybrid_evolver.hpp"
#include "labelled/graphics_interface.hpp"
#include "labelled/figure.hpp"
#include "hybrid/hybrid_graphics.hpp"
#include "analysis/hybrid_reachability_analyser.hpp"
#include "drawing/drawer_manager.hpp"
#include "drawing/drawer.hpp"
#include "drawing/drawing_command_line_interface.hpp"

#include "macro/testing.hpp"

using namespace Ariadne;
using namespace std;


Colour reach_set_colour(0.25,0.25,0.50);
Colour intermediate_set_colour(0.50,0.50,0.75);
Colour final_set_colour(0.75,0.75,1.00);
Colour initial_set_colour(0.75,0.75,1.00);
Colour guard_set_colour(0.75,0.75,0.75);

// Test evolution of realistic hybrid systems
class TestHybridReachability
{
    mutable shared_ptr<HybridReachabilityAnalyser> analyser;
  private:
    Void _set_analyser(const HybridAutomatonInterface& system) const;
  public:
    Void test() const;
    Void test_bouncing_ball() const;
};

Void TestHybridReachability::_set_analyser(const HybridAutomatonInterface& system) const
{
    GeneralHybridEvolver evolver(system);
    evolver.set_integrator(TaylorPicardIntegrator(StepMaximumError(1e-5)));
    evolver.configuration().set_maximum_step_size(0.25);
    evolver.configuration().set_maximum_enclosure_radius(0.125);
    evolver.configuration().set_maximum_enclosure_radius(0.5);
    analyser.reset(new HybridReachabilityAnalyser(evolver));
    analyser->configuration().set_maximum_grid_fineness(2);
    analyser->configuration().set_lock_to_grid_time(5);
    ARIADNE_TEST_PRINT(analyser->configuration());
}

Void TestHybridReachability::test() const {
    ARIADNE_TEST_CALL(test_bouncing_ball());
}

Void TestHybridReachability::test_bouncing_ball() const {
    HybridAutomaton bouncing_ball;
    Real one(1);
    DiscreteLocation q;
    DiscreteEvent e("e");
    RealVariable x("x");
    RealVariable v("v");
    TimeVariable t;

    Real lambda(0.5_x);
    bouncing_ball.new_mode(q,{dot(x)=v,dot(v)=-one});
    bouncing_ball.new_transition(q,e,q,{next(x)=x,next(v)=-lambda*v},x<=0,EventKind::IMPACT);
    ARIADNE_TEST_PRINT(bouncing_ball);

    Dyadic height=2;
    Dyadic radius=1/pow(two,6);
    HybridBoundedConstraintSet initial(q,{height-radius<=x<=height+radius,-radius<=v<=radius});
    // Exercise the first impact, not a long sequence of bounces.
    Decimal tmax(2.5);
    Natural maxsteps=2u;
    HybridTime time(tmax,maxsteps);

    this->_set_analyser(bouncing_ball);

    auto reach_evolve=analyser->upper_reach_evolve(initial,time);

    ARIADNE_TEST_ASSERT(reach_evolve.first.size()>0u);
    ARIADNE_TEST_ASSERT(reach_evolve.second.size()>0u);
}


Int main(Int argc, const char* argv[])
{
    if (not DrawingCommandLineInterface().acquire(argc,argv)) return -1;

    DrawerManager::instance().set_drawer(GridDrawer(2));

    TestHybridReachability().test();
    return ARIADNE_TEST_FAILURES;
}

