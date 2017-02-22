/***************************************************************************
 *            test_hybrid_reachability_analysis.cpp
 *
 *  Copyright  2006-8  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <iostream>
#include <fstream>
#include <string>

#include "config.h"
#include "function/taylor_function.hpp"
#include "function/formula.hpp"
#include "algebra/algebra.hpp"
#include "geometry/function_set.hpp"
#include "geometry/grid_set.hpp"
#include "hybrid/hybrid_time.hpp"
#include "hybrid/hybrid_set.hpp"
#include "hybrid/hybrid_automata.hpp"
#include "hybrid/hybrid_evolver.hpp"
#include "hybrid/hybrid_reachability_analyser.hpp"
#include "hybrid/hybrid_graphics.hpp"
#include "utility/logging.hpp"

#include "test.hpp"

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


namespace Ariadne {

HybridExactBoxes
bounding_boxes(const HybridSpaceInterface& space, ExactIntervalType bound)
{
    HybridExactBoxes result;
    Set<DiscreteLocation> locations = dynamic_cast<const MonolithicHybridSpace&>(space).locations();
    for(Set<DiscreteLocation>::ConstIterator loc_iter=locations.begin();
        loc_iter!=locations.end(); ++loc_iter)
        {
            DiscreteLocation const& loc=*loc_iter;
            result.insert(loc,space[loc],ExactBoxType(space[loc].dimension(), bound));
        }
    return result;
}

}


class TestHybridReachabilityAnalyser
{

    typedef GeneralHybridEvolver HybridEvolverType;
    typedef HybridEvolverType::EnclosureType HybridEnclosureType;
    typedef HybridEnclosureType::ContinuousStateSetType ContinuousEnclosureType;

    shared_ptr<HybridReachabilityAnalyser> analyser_ptr;
    HybridAutomaton system;
    HybridReachabilityAnalyser analyser;
    Grid grid;
    ExactIntervalType bound;
    HybridSet initial_set;
    HybridSet safe_set;
    HybridTime reach_time;
    Axes2d axes;
  public:

    static HybridAutomaton build_system()
    {
        HybridAutomaton sys;
        cout << "Done initialising variables\n";
        std::cout<<std::setprecision(20);
        std::cerr<<std::setprecision(20);
        std::clog<<std::setprecision(20);
        DiscreteLocation location(1);
        RealVariable x("x");
        RealVariable y("y");
        Real a(-0.25); Real b(1.0);
        sys.new_mode(location,{dot(x)=-a*x-b*y,dot(y)=b*x+2*a*y} );
        cout << "Done building system\n";
        return sys;
    }

    static HybridReachabilityAnalyser build_analyser(const HybridAutomatonInterface& system)
    {
        TaylorFunctionFactory function_factory(ThresholdSweeper<Float64>(Precision64(),1e-8));
        GeneralHybridEvolver evolver(system,function_factory);

        HybridReachabilityAnalyser analyser(system,evolver);
        analyser.configuration().set_maximum_grid_depth(4);
        analyser.configuration().set_lock_to_grid_time(1.0);
        cout << "Done building analyser\n";
        return analyser;
    }

    TestHybridReachabilityAnalyser(Nat analyser_verbosity = 0u)
        : system(build_system()),
          analyser(build_analyser(system)),
          grid(2),
          bound(-4,4),
          reach_time(3.0,4),
          axes({-3<=RealVariable("x")<=3,-3<=RealVariable("y")<=3})
    {
        analyser.verbosity=analyser_verbosity;
        DiscreteLocation location(1);

        RealVariable x("x");
        RealVariable y("y");

        //ImageSet initial_box(make_box("[1.99,2.01]x[-0.01,0.01]"));
        Decimal x0l(2.01), x0u(2.02), y0l(0.01), y0u(0.02);
        initial_set=HybridSet(location,{x.in(x0l,(Real)x0u),y.in(y0l,(Real)y0u)});
        safe_set=HybridSet(location,{x.in(-2,(Real)3),y.in(-2,(Real)3)},{sqr(x)+sqr(y)<=sqr((Real)3)});
        cout << "Done creating initial and safe sets\n" << endl;

        cout << "system=" << system << endl;
        cout << "initial_set=" << initial_set << endl;
        cout << "safe_set=" << safe_set << endl;

        //ARIADNE_ASSERT(inside(initial_set[loc],bounding_set[loc]));

    }

    Void test_lower_reach_lower_evolve() {
        DiscreteLocation loc(1);
        ExactBoxType bounding_box(2,bound);
        cout << "Computing timed reachable set" << endl;
        HybridGridTreeSet hybrid_lower_reach=analyser.lower_reach(initial_set,reach_time);
        GridTreeSet& lower_reach=hybrid_lower_reach[loc];

        ARIADNE_TEST_ASSERT(lower_reach.size() > 0);

        cout << "Computing timed evolve set" << endl;
        HybridGridTreeSet hybrid_lower_evolve=analyser.lower_evolve(initial_set,reach_time);
        GridTreeSet& lower_evolve=hybrid_lower_evolve[loc];

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
        Pair<HybridGridTreeSet,HybridGridTreeSet> reach_evolve_set = analyser.lower_reach_evolve(initial_set,reach_time);
        GridTreeSet& lower_reach=reach_evolve_set.first[loc];
        GridTreeSet& lower_evolve=reach_evolve_set.second[loc];
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
        HybridGridTreeSet upper_reach_set=analyser.upper_reach(initial_set,reach_time);
        cout << "upper_reach_set="<<upper_reach_set<<std::endl;

        ARIADNE_TEST_ASSERT(upper_reach_set.size() > 0);

        cout << "Computing timed evolve set" << endl;
        HybridGridTreeSet upper_evolve_set=analyser.upper_evolve(initial_set,reach_time);
        cout << "upper_evolve_set="<<upper_evolve_set<<std::endl;

        ARIADNE_TEST_ASSERT(upper_evolve_set.size() > 0);

        plot("test_reachability_analyser-map_upper_reach_upper_evolve.png",axes,
             reach_set_colour,upper_reach_set,evolve_set_colour,upper_evolve_set);
    }

    Void test_upper_reach_evolve() {
        cout << "Computing timed reach-evolve set" << endl;
        DiscreteLocation loc(1);
        ExactBoxType bounding_box(2,bound);
        Pair<HybridGridTreeSet,HybridGridTreeSet> reach_evolve_set = analyser.upper_reach_evolve(initial_set,reach_time);

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
        HybridGridTreeSet lower_reach_set=analyser.lower_reach(initial_set);

        ARIADNE_TEST_ASSERT(lower_reach_set.size() > 0);

        plot("test_reachability_analyser-map_infinite_time_lower_reach.png",axes,
             reach_set_colour,lower_reach_set);
    }

    Void test_outer_chain_reach() {
        cout << "Computing outer chain reachable set" << endl;
        DiscreteLocation loc(1);
        HybridExactBoxes bounding_boxes
            =Ariadne::bounding_boxes(system.state_space(),bound);

        analyser.configuration().set_transient_steps(1);
        analyser.configuration().set_transient_time(12.0);
        analyser.configuration().set_lock_to_grid_time(6.0);
        analyser.configuration().set_maximum_grid_depth(3);
        analyser.configuration().set_bounding_domain_ptr(shared_ptr<HybridExactBoxes>(new HybridExactBoxes(bounding_boxes)));
        cout << analyser.configuration();

        HybridGridTreeSet outer_chain_reach_set=analyser.outer_chain_reach(initial_set);

        ARIADNE_TEST_ASSERT(outer_chain_reach_set.size() > 0);

        plot("test_reachability_analyser-map_outer_chain_reach.png",axes,
             reach_set_colour,outer_chain_reach_set);

        cout << "Recomputing with tight restriction" << endl;

        analyser.configuration().set_maximum_grid_height(1);
        ARIADNE_TEST_THROWS(analyser.outer_chain_reach(initial_set),OuterChainOverspill);
    }

    Void test_verify_safety() {
        cout << "Verifying safety" << endl;
        DiscreteLocation loc(1);

        cout << analyser.configuration();

        auto safety_certificate=analyser.verify_safety(initial_set,safe_set);

        ARIADNE_TEST_ASSERT(definitely(safety_certificate.is_safe));

        plot("test_reachability_analyser-verify_safety.png",axes,
             safe_set_colour,safe_set,
             reach_set_colour,safety_certificate.chain_reach_set,
             initial_set_colour,initial_set);

    }

    Void test() {
        //ValidatedTaylorModel64::set_default_sweep_threshold(1e-6);
        //ValidatedTaylorModel64::set_default_maximum_degree(6u);

        ARIADNE_TEST_CALL(test_lower_reach_lower_evolve());
        ARIADNE_TEST_CALL(test_lower_reach_evolve());
        ARIADNE_TEST_CALL(test_upper_reach_upper_evolve());
        ARIADNE_TEST_CALL(test_upper_reach_evolve());
        ARIADNE_TEST_CALL(test_infinite_time_lower_reach());

        ARIADNE_TEST_CALL(test_outer_chain_reach());
        ARIADNE_TEST_CALL(test_verify_safety());
    }

};


Int main(Int argc, const char* argv[])
{
    Int analyser_verbosity=get_verbosity(argc,argv);

    TestHybridReachabilityAnalyser(analyser_verbosity).test();
    return ARIADNE_TEST_FAILURES;
}

