/***************************************************************************
 *            test_reachability_analysis.cpp
 *
 *  Copyright  2006-8  Pieter Collins
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

#include <iostream>
#include <fstream>
#include <string>

#include "config.hpp"
#include "function/taylor_function.hpp"
#include "function/formula.hpp"
#include "algebra/algebra.hpp"
#include "geometry/function_set.hpp"
#include "geometry/grid_paving.hpp"
#include "dynamics/vector_field_evolver.hpp"
#include "dynamics/reachability_analyser.hpp"
#include "symbolic/expression_set.hpp"
#include "solvers/integrator.hpp"
#include "output/colour.hpp"
#include "output/logging.hpp"

#include "../test.hpp"

using namespace std;
using namespace Ariadne;

Colour reach_set_colour(0.25,0.25,0.50);
Colour intermediate_set_colour(0.50,0.50,0.75);
Colour final_set_colour(0.75,0.75,1.00);
Colour evolve_set_colour(0.75,0.75,1.00);
Colour initial_set_colour(0.75,0.75,1.00);
Colour guard_set_colour(0.75,0.75,0.75);
Colour bounding_box_colour(0.75,0.75,0.875);
Colour safe_set_colour(0.50,1.00,0.50);

class TestReachabilityAnalyser
{

    typedef VectorFieldEvolver EvolverType;
    typedef EvolverType::SystemType SystemType;
    typedef SystemType::StateSpaceType StateSpaceType;
    typedef EvolverType::EnclosureType EnclosureType;
    typedef RealExpressionBoundedConstraintSet SetType;
    typedef ReachabilityAnalyser<SystemType> AnalyserType;

    SystemType system;
    AnalyserType analyser;
    Grid grid;
    RealVariable x;
    RealVariable y;
    ExactBoxType bounding;
    SetType initial_set;
    SetType safe_set;
    Real reach_time;
  public:

    static SystemType build_system()
    {
        std::cout<<std::setprecision(20);
        std::cerr<<std::setprecision(20);
        std::clog<<std::setprecision(20);
        RealVariable x("x");
        RealVariable y("y");
        Real a(-0.5); Real b(1.0);
        SystemType sys({dot(x)=-a*x-b*y,dot(y)=b*x+2*a*y});
        cout << "Done building system\n";
        return sys;
    }

    static AnalyserType build_analyser(const SystemType& system)
    {
        TaylorPicardIntegrator integrator(maximum_error=1e-4,sweep_threshold=1e-8,lipschitz_constant=0.5,
                                                 step_maximum_error=1e-6,step_sweep_threshold=1e-10,maximum_temporal_order=8);
        EvolverType evolver(system,integrator);

        AnalyserType analyser(system,evolver);
        analyser.configuration().set_maximum_grid_depth(4);
        cout << "Done building analyser\n";
        return analyser;
    }

    TestReachabilityAnalyser(Nat analyser_verbosity = 0u)
        : system(build_system()),
          analyser(build_analyser(system)),
          grid(2),
          x("x"),
          y("y"),
          bounding({{-4,4},{-4,4}}),
          initial_set({1.98_dec<=x<=1.99_dec,0.01_dec<=y<=0.02_dec}),
          safe_set({-2_dec<=x<=3_dec,-2_dec<=y<=3_dec}),//,{sqr(x)+sqr(y)<=sqr((Real)3)}),
          reach_time(3.0)
    {
        analyser.verbosity=analyser_verbosity;

        cout << "Done creating initial and safe sets\n" << endl;

        cout << "system=" << system << endl;
        cout << "initial_set=" << initial_set << endl;
        cout << "safe_set=" << safe_set << endl;
    }
/*
    Void test_lower_reach_lower_evolve() {
        DiscreteLocation loc(1);
        ExactBoxType bounding_box(2,bound);
        cout << "Computing timed reachable set" << endl;
        HybridGridTreePaving hybrid_lower_reach=analyser.lower_reach(initial_set,reach_time);
        GridTreePaving& lower_reach=hybrid_lower_reach[loc];

        ARIADNE_TEST_ASSERT(lower_reach.size() > 0);

        cout << "Computing timed evolve set" << endl;
        HybridGridTreePaving hybrid_lower_evolve=analyser.lower_evolve(initial_set,reach_time);
        GridTreePaving& lower_evolve=hybrid_lower_evolve[loc];

        ARIADNE_TEST_ASSERT(lower_evolve.size() > 0);

        cout << "Reached " << lower_reach.size() << " cells " << endl << endl;
        cout << "Evolved to " << lower_evolve.size() << " cells " << endl << endl;

        plot("test_reachability_analyser-map_lower_reach_lower_evolve.png",axes,
             reach_set_colour,hybrid_lower_reach,evolve_set_colour,hybrid_lower_evolve);
    }

    Void test_lower_reach_evolve() {
        DiscreteLocation loc(1);
        ExactBoxType bounding_box(2,bound);

        cout << "Computing timed reach-evolve set" << endl;
        Pair<HybridGridTreePaving,HybridGridTreePaving> reach_evolve_set = analyser.lower_reach_evolve(initial_set,reach_time);
        GridTreePaving& lower_reach=reach_evolve_set.first[loc];
        GridTreePaving& lower_evolve=reach_evolve_set.second[loc];
        cout << "Reached " << lower_reach.size() << " cells " << endl << endl;
        cout << "Evolved to " << lower_evolve.size() << " cells " << endl << endl;

        ARIADNE_TEST_ASSERT(lower_reach.size() > 0);
        ARIADNE_TEST_ASSERT(lower_evolve.size() > 0);

        plot("test_reachability_analyser-map_lower_reach_evolve.png",axes,
             reach_set_colour,reach_evolve_set.first,evolve_set_colour,reach_evolve_set.second);
    }

    Void test_upper_reach_upper_evolve() {
        DiscreteLocation loc(1);
        ExactBoxType bounding_box(2,bound);
        cout << "Computing timed reachable set" << endl;
        HybridGridTreePaving upper_reach_set=analyser.upper_reach(initial_set,reach_time);
        cout << "upper_reach_set="<<upper_reach_set<<std::endl;

        ARIADNE_TEST_ASSERT(upper_reach_set.size() > 0);

        cout << "Computing timed evolve set" << endl;
        HybridGridTreePaving upper_evolve_set=analyser.upper_evolve(initial_set,reach_time);
        cout << "upper_evolve_set="<<upper_evolve_set<<std::endl;

        ARIADNE_TEST_ASSERT(upper_evolve_set.size() > 0);

        plot("test_reachability_analyser-map_upper_reach_upper_evolve.png",axes,
             reach_set_colour,upper_reach_set,evolve_set_colour,upper_evolve_set);
    }

    Void test_upper_reach_evolve() {
        cout << "Computing timed reach-evolve set" << endl;
        DiscreteLocation loc(1);
        ExactBoxType bounding_box(2,bound);
        Pair<HybridGridTreePaving,HybridGridTreePaving> reach_evolve_set = analyser.upper_reach_evolve(initial_set,reach_time);

        ARIADNE_TEST_ASSERT(reach_evolve_set.first.size() > 0);
        ARIADNE_TEST_ASSERT(reach_evolve_set.second.size() > 0);

        plot("test_reachability_analyser-map_upper_reach_evolve.png",axes,
             reach_set_colour,reach_evolve_set.first,final_set_colour,reach_evolve_set.second);
    }

    Void test_infinite_time_lower_reach() {

        DiscreteLocation loc(1);
        HybridExactBoxes bounding_boxes
            =Ariadne::bounding_boxes(system.state_space(),bound);

        analyser.configuration().set_transient_time(4.0);
        analyser.configuration().set_bounding_domain_ptr(shared_ptr<HybridExactBoxes>(new HybridExactBoxes(bounding_boxes)));
        cout << analyser.configuration();

        cout << "Computing infinite time lower reachable set" << endl;
        HybridGridTreePaving lower_reach_set=analyser.lower_reach(initial_set);

        ARIADNE_TEST_ASSERT(lower_reach_set.size() > 0);

        plot("test_reachability_analyser-map_infinite_time_lower_reach.png",axes,
             reach_set_colour,lower_reach_set);
    }

    Void test_outer_chain_reach() {
        cout << "Computing outer chain reachable set" << endl;

        analyser.configuration().set_transient_steps(1);
        analyser.configuration().set_transient_time(12.0);
        analyser.configuration().set_lock_to_grid_time(6.0);
        analyser.configuration().set_maximum_grid_depth(3);
        analyser.configuration().set_bounding_domain_ptr(shared_ptr<ExactBox>(new ExactBox(system.,bound)));
        cout << analyser.configuration();

        HybridGridTreePaving outer_chain_reach_set=analyser.outer_chain_reach(initial_set);

        ARIADNE_TEST_ASSERT(outer_chain_reach_set.size() > 0);

        plot("test_reachability_analyser-map_outer_chain_reach.png",axes,
             reach_set_colour,outer_chain_reach_set);

        cout << "Recomputing with tight restriction" << endl;

        analyser.configuration().set_maximum_grid_height(1);
        ARIADNE_TEST_THROWS(analyser.outer_chain_reach(initial_set),OuterChainOverspill);
    }

    Void test_verify_safety() {
        cout << "Verifying safety" << endl;

        cout << analyser.configuration();

        auto safety_certificate=analyser.verify_safety(initial_set,safe_set);

        ARIADNE_TEST_ASSERT(definitely(safety_certificate.is_safe));

        auto grid=analyser.configuration().grid();
        auto safe_cells=inner_approximation(safe_set, grid, analyser.configuration().maximum_grid_depth());
        plot("test_reachability_analyser-verify_safety.png",axes,
             safe_set_colour,safety_certificate.safe_set,
             reach_set_colour,safety_certificate.chain_reach_set,
             initial_set_colour,initial_set);

    }
*/
    Void test() {
/*
        ARIADNE_TEST_CALL(test_lower_reach_lower_evolve());
        ARIADNE_TEST_CALL(test_lower_reach_evolve());
        ARIADNE_TEST_CALL(test_upper_reach_upper_evolve());
        ARIADNE_TEST_CALL(test_upper_reach_evolve());
        ARIADNE_TEST_CALL(test_infinite_time_lower_reach());

        ARIADNE_TEST_CALL(test_outer_chain_reach());
        ARIADNE_TEST_CALL(test_verify_safety());
*/
    }

};


Int main(Int argc, const char* argv[])
{
    unsigned int analyser_verbosity=get_verbosity(argc,argv);

    TestReachabilityAnalyser(analyser_verbosity).test();
    return ARIADNE_TEST_FAILURES;
}

