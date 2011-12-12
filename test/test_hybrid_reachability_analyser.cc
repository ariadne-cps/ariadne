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
#include "evolution_parameters.h"
#include "hybrid_evolver.h"
#include "hybrid_discretiser.h"
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

    HybridReachabilityAnalyser analyser;
    HybridAutomaton system;
    Grid grid;
    Interval bound;
    HybridSet initial_set;
    HybridTime reach_time;
	
  public:
    static HybridReachabilityAnalyser build_analyser()
    {
        EvolutionParameters parameters;
        parameters.maximum_grid_depth=4;
        parameters.maximum_step_size=0.25;
        parameters.lock_to_grid_time=1.0;
        TaylorFunctionFactory factory(ThresholdSweeper(1e-8));

        Grid grid(2);
        HybridEvolverType evolver(parameters,factory);
        //HybridDiscretiser<EnclosureType> discretiser(evolver);
        HybridReachabilityAnalyser analyser(parameters,evolver);
        cout << "Done building analyser\n";
        return analyser;
    }

    TestHybridReachabilityAnalyser()
        : analyser(build_analyser()),
          system(),
          grid(2),
          bound(-4,4),
          reach_time(3.0,4)
    {
        cout << "Done initialising variables\n";
        std::cout<<std::setprecision(20);
        std::cerr<<std::setprecision(20);
        std::clog<<std::setprecision(20);
        DiscreteLocation location(StringVariable("q")|"1");

        RealVariable x("x");
        RealVariable y("y");
        system.new_mode(location,(dot(x)=-0.5*x-1.0*y,dot(y)=1.0*x-0.5*y) );
        cout << "Done building system\n";

        //ImageSet initial_box(make_box("[1.99,2.01]x[-0.01,0.01]"));
        initial_set=HybridSet(location,(x.in(2.01,2.02),y.in(0.01,0.02)));
        cout << "Done creating initial set\n" << endl;

        cout << "system=" << system << endl;
        cout << "initial_set=" << initial_set << endl;

        //ARIADNE_ASSERT(inside(initial_set[loc],bounding_set[loc]));

    }

    void test_lower_reach_evolve() {
        DiscreteLocation loc(1);
        Box bounding_box(2,bound);
        cout << "Computing timed evolve set" << endl;
        HybridGridTreeSet hybrid_lower_evolve=analyser.lower_evolve(system,initial_set,reach_time);
        cout << "Computing timed reachable set" << endl;
        HybridGridTreeSet hybrid_lower_reach=analyser.lower_reach(system,initial_set,reach_time);
        GridTreeSet& lower_evolve=hybrid_lower_evolve[loc];
        GridTreeSet& lower_reach=hybrid_lower_reach[loc];
        RealBoundedConstraintSet const& initial=initial_set[loc];
        cout << "Evolved to " << lower_evolve.size() << " cells " << endl << endl;
        cout << "Reached " << lower_reach.size() << " cells " << endl << endl;
        RealVariable x("x"); RealVariable y("y");
   		Axes2d axes(-1.0<=x<=+3.0,-2.0<=y<=+2.0);
        plot("test_reachability_analyser-map_lower_reach_evolve.png",axes,
             reach_set_colour,hybrid_lower_reach,evolve_set_colour,hybrid_lower_evolve);
    }

    void test_upper_reach_evolve() {
        cout << "Computing timed reachable set" << endl;
        DiscreteLocation loc(1);
        Box bounding_box(2,bound);
        HybridGridTreeSet upper_evolve_set=analyser.upper_evolve(system,initial_set,reach_time);
        cout << "upper_evolve_set="<<upper_evolve_set<<std::endl;
        HybridGridTreeSet upper_reach_set=analyser.upper_reach(system,initial_set,reach_time);
        cout << "upper_reach_set="<<upper_reach_set<<std::endl;

        const GridTreeSet& upper_evolve=upper_evolve_set[loc];
        const GridTreeSet& upper_reach=upper_reach_set[loc];
        RealBoundedConstraintSet const& initial=initial_set[loc];
        RealVariable x("x"); RealVariable y("y");
   		Axes2d axes(-1.0<=x<=+3.0,-2.0<=y<=+2.0);
        plot("test_reachability_analyser-map_upper_reach_evolve.png",axes,
             reach_set_colour,upper_reach_set,final_set_colour,upper_evolve_set);
    }

    void test_chain_reach() {
        cout << "Computing chain reachable set" << endl;
        DiscreteLocation loc(1);
        HybridBoxes bounding_boxes
            =Ariadne::bounding_boxes(system.state_space(),bound);
        Box bounding_box=bounding_boxes[loc];

        analyser.verbosity=0;
        analyser.parameters().transient_time=4.0;
        cout << analyser.parameters();
        HybridGridTreeSet chain_reach_set=analyser.chain_reach(system,initial_set,bounding_boxes);
        RealVariable x("x"); RealVariable y("y");
   		Axes2d axes(-1.0<=x<=+3.0,-2.0<=y<=+2.0);
        plot("test_reachability_analyser-map_chain_reach.png",axes,
             reach_set_colour,chain_reach_set);
    }

    void test() {
        //IntervalTaylorModel::set_default_sweep_threshold(1e-6);
        //IntervalTaylorModel::set_default_maximum_degree(6u);

        ARIADNE_TEST_CALL(test_lower_reach_evolve());
        ARIADNE_TEST_CALL(test_upper_reach_evolve());
        ARIADNE_TEST_CALL(test_chain_reach());
    }

};


int main(int nargs, const char* args[])
{
    TestHybridReachabilityAnalyser().test();
    return ARIADNE_TEST_FAILURES;
}

