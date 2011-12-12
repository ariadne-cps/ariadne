/***************************************************************************
 *            test_hybrid_reachability_analysis.cc
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
#include "taylor_function.h"
#include "function_set.h"
#include "grid_set.h"
#include "hybrid_time.h"
#include "hybrid_set.h"
#include "hybrid_automaton.h"
#include "hybrid_evolver.h"
#include "hybrid_reachability_analyser.h"
#include "hybrid_graphics.h"
#include "logging.h"

#include "test.h"

using namespace std;
using namespace Ariadne;

Colour reach_set_colour(0.25,0.25,0.50);
Colour intermediate_set_colour(0.50,0.50,0.75);
Colour final_set_colour(0.75,0.75,1.00);
Colour evolve_set_colour(0.75,0.75,1.00);
Colour initial_set_colour(0.75,0.75,1.00);
Colour guard_set_colour(0.75,0.75,0.75);
Colour bounding_box_colour(0.75,0.75,0.875);


namespace Ariadne {

HybridBoxes
bounding_boxes(const HybridSpaceInterface& space, Interval bound)
{
    HybridBoxes result;
    Set<DiscreteLocation> locations = dynamic_cast<const MonolithicHybridSpace&>(space).locations();
    for(Set<DiscreteLocation>::const_iterator loc_iter=locations.begin();
        loc_iter!=locations.end(); ++loc_iter)
        {
            result.insert(make_pair(*loc_iter,Box(space[*loc_iter].dimension(), bound)));
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
    Interval bound;
    HybridSet initial_set;
    HybridTime reach_time;
	
  public:

    static HybridAutomaton build_system()
    {
        HybridAutomaton sys;
        cout << "Done initialising variables\n";
        std::cout<<std::setprecision(20);
        std::cerr<<std::setprecision(20);
        std::clog<<std::setprecision(20);
        DiscreteLocation location(StringVariable("q")|"1");

        RealVariable x("x");
        RealVariable y("y");
        sys.new_mode(location,(dot(x)=-0.5*x-1.0*y,dot(y)=1.0*x-0.5*y) );
        cout << "Done building system\n";
        return sys;
    }

    static HybridReachabilityAnalyser build_analyser(const HybridAutomatonInterface& system)
    {
        TaylorFunctionFactory function_factory(ThresholdSweeper(1e-8));
        GeneralHybridEvolverFactory evolver_factory(function_factory);

        HybridReachabilityAnalyser analyser(system,evolver_factory);
        analyser.configuration().set_maximum_grid_depth(4);
        analyser.configuration().set_lock_to_grid_time(1.0);
        cout << "Done building analyser\n";
        return analyser;
    }

    TestHybridReachabilityAnalyser()
        : system(build_system()),
          analyser(build_analyser(system)),
          grid(2),
          bound(-4,4),
          reach_time(3.0,4)
    {
        DiscreteLocation location(StringVariable("q")|"1");

        RealVariable x("x");
        RealVariable y("y");

        //ImageSet initial_box(make_box("[1.99,2.01]x[-0.01,0.01]"));
        initial_set=HybridSet(location,(x.in(2.01,2.02),y.in(0.01,0.02)));
        cout << "Done creating initial set\n" << endl;

        cout << "system=" << system << endl;
        cout << "initial_set=" << initial_set << endl;

        //ARIADNE_ASSERT(inside(initial_set[loc],bounding_set[loc]));

    }

    void test_lower_reach_lower_evolve() {
        DiscreteLocation loc(1);
        Box bounding_box(2,bound);
        cout << "Computing timed reachable set" << endl;
        HybridGridTreeSet hybrid_lower_reach=analyser.lower_reach(initial_set,reach_time);
        GridTreeSet& lower_reach=hybrid_lower_reach[loc];
        cout << "Computing timed evolve set" << endl;
        HybridGridTreeSet hybrid_lower_evolve=analyser.lower_evolve(initial_set,reach_time);
        GridTreeSet& lower_evolve=hybrid_lower_evolve[loc];
        cout << "Reached " << lower_reach.size() << " cells " << endl << endl;
        cout << "Evolved to " << lower_evolve.size() << " cells " << endl << endl;
        RealVariable x("x"); RealVariable y("y");
        Axes2d axes(-1.0<=x<=+3.0,-2.0<=y<=+2.0);
        plot("test_reachability_analyser-map_lower_reach_lower_evolve.png",axes,
             reach_set_colour,hybrid_lower_reach,evolve_set_colour,hybrid_lower_evolve);
    }

    void test_lower_reach_evolve() {
        DiscreteLocation loc(1);
        Box bounding_box(2,bound);

        cout << "Computing timed reach-evolve set" << endl;
        Pair<HybridGridTreeSet,HybridGridTreeSet> result = analyser.lower_reach_evolve(initial_set,reach_time);
        GridTreeSet& lower_reach=result.first[loc];
        GridTreeSet& lower_evolve=result.second[loc];
        cout << "Reached " << lower_reach.size() << " cells " << endl << endl;
        cout << "Evolved to " << lower_evolve.size() << " cells " << endl << endl;
        RealVariable x("x"); RealVariable y("y");
        Axes2d axes(-1.0<=x<=+3.0,-2.0<=y<=+2.0);
        plot("test_reachability_analyser-map_lower_reach_evolve.png",axes,
             reach_set_colour,result.first,evolve_set_colour,result.second);
    }

    void test_upper_reach_upper_evolve() {
        DiscreteLocation loc(1);
        Box bounding_box(2,bound);
        cout << "Computing timed reachable set" << endl;
        HybridGridTreeSet upper_reach_set=analyser.upper_reach(initial_set,reach_time);
        cout << "upper_reach_set="<<upper_reach_set<<std::endl;
        cout << "Computing timed evolve set" << endl;
        HybridGridTreeSet upper_evolve_set=analyser.upper_evolve(initial_set,reach_time);
        cout << "upper_evolve_set="<<upper_evolve_set<<std::endl;

        RealVariable x("x"); RealVariable y("y");
        Axes2d axes(-1.0<=x<=+3.0,-2.0<=y<=+2.0);
        plot("test_reachability_analyser-map_upper_reach_upper_evolve.png",axes,
             reach_set_colour,upper_reach_set,evolve_set_colour,upper_evolve_set);
    }

    void test_upper_reach_evolve() {
        cout << "Computing timed reach-evolve set" << endl;
        DiscreteLocation loc(1);
        Box bounding_box(2,bound);
        Pair<HybridGridTreeSet,HybridGridTreeSet> result = analyser.upper_reach_evolve(initial_set,reach_time);

        RealVariable x("x"); RealVariable y("y");
        Axes2d axes(-1.0<=x<=+3.0,-2.0<=y<=+2.0);
        plot("test_reachability_analyser-map_upper_reach_evolve.png",axes,
             reach_set_colour,result.first,final_set_colour,result.second);
    }

    void test_infinite_time_lower_reach() {

        DiscreteLocation loc(1);
        HybridBoxes bounding_boxes
            =Ariadne::bounding_boxes(system.state_space(),bound);
        Box bounding_box=bounding_boxes[loc];

        analyser.configuration().set_transient_time(4.0);
        analyser.configuration().set_bounding_domain_ptr(shared_ptr<HybridBoxes>(new HybridBoxes(bounding_boxes)));
        cout << analyser.configuration();

        cout << "Computing infinite time lower reachable set" << endl;

        HybridGridTreeSet lower_reach_set=analyser.lower_reach(initial_set);
        RealVariable x("x"); RealVariable y("y");
        Axes2d axes(-1.0<=x<=+3.0,-2.0<=y<=+2.0);
        plot("test_reachability_analyser-map_infinite_time_lower_reach.png",axes,
             reach_set_colour,lower_reach_set);
    }

    void test_outer_chain_reach() {
        cout << "Computing outer chain reachable set" << endl;
        DiscreteLocation loc(1);
        HybridBoxes bounding_boxes
            =Ariadne::bounding_boxes(system.state_space(),bound);
        Box bounding_box=bounding_boxes[loc];

        analyser.configuration().set_transient_time(4.0);
        analyser.configuration().set_bounding_domain_ptr(shared_ptr<HybridBoxes>(new HybridBoxes(bounding_boxes)));
        cout << analyser.configuration();

        HybridGridTreeSet outer_chain_reach_set=analyser.outer_chain_reach(initial_set);
        RealVariable x("x"); RealVariable y("y");
        Axes2d axes(-1.0<=x<=+3.0,-2.0<=y<=+2.0);
        plot("test_reachability_analyser-map_outer_chain_reach.png",axes,
             reach_set_colour,outer_chain_reach_set);

        cout << "Recomputing with tight restriction" << endl;

        analyser.configuration().set_maximum_grid_height(1);
        ARIADNE_TEST_THROWS(analyser.outer_chain_reach(initial_set),OuterChainOverspill);


    }

    void test() {
        //IntervalTaylorModel::set_default_sweep_threshold(1e-6);
        //IntervalTaylorModel::set_default_maximum_degree(6u);

        ARIADNE_TEST_CALL(test_lower_reach_lower_evolve());
        ARIADNE_TEST_CALL(test_lower_reach_evolve());
        ARIADNE_TEST_CALL(test_upper_reach_upper_evolve());
        ARIADNE_TEST_CALL(test_upper_reach_evolve());
        ARIADNE_TEST_CALL(test_infinite_time_lower_reach());
        ARIADNE_TEST_CALL(test_outer_chain_reach());
    }

};


int main(int nargs, const char* args[])
{
    TestHybridReachabilityAnalyser().test();
    return ARIADNE_TEST_FAILURES;
}

