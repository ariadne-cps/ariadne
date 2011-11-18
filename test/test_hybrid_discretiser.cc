/***************************************************************************
 *            test_hybrid_discretiser.cc
 *
 *  Copyright  2006-11  Pieter Collins
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

#include <fstream>

#include "config.h"
#include "tuple.h"
#include "vector.h"
#include "matrix.h"
#include "function.h"
#include "taylor_set.h"
#include "zonotope.h"
#include "list_set.h"
#include "grid_set.h"
#include "integrator.h"
#include "hybrid_evolver.h"
#include "hybrid_time.h"
#include "hybrid_space.h"
#include "hybrid_set.h"
#include "hybrid_automaton.h"
#include "hybrid_orbit.h"
#include "hybrid_discretiser.h"

#include "graphics.h"
#include "logging.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

class TestHybridDiscretiser
{
  public:
    void test() const;
  private:
    void test_hybrid_time() const;
};

int main()
{
    TestHybridDiscretiser().test();
    return ARIADNE_TEST_FAILURES;
}

void TestHybridDiscretiser::test() const
{
    ARIADNE_TEST_CALL(test_hybrid_time());
}



void TestHybridDiscretiser::test_hybrid_time() const
{
    typedef GeneralHybridEvolver EvolverType;
    typedef EvolverType::EnclosureType EnclosureType;
    typedef EnclosureType::ContinuousStateSetType ContinuousEnclosureType;

    // Set up the evolution parameters and grid
    Real time(2.5);
    uint steps(6);
    double maximum_step_size(4.0);
    int depth=3;
    DiscreteLocation location1(1);
    DiscreteLocation location2(2);
    DiscreteLocation location3(3);
    DiscreteLocation location4(4);
    DiscreteEvent event1(1);
    DiscreteEvent event2(2);

    EvolutionParameters parameters;
    parameters.maximum_step_size=maximum_step_size;
    TaylorFunctionFactory factory(ThresholdSweeper(1e-8));
    Grid grid(2);

    // Set up the evaluators
    EvolverType evolver(parameters,factory);
    HybridDiscretiser< EnclosureType > discrete_evolver(evolver);
    discrete_evolver.verbosity=9;

    // Set up the vector field
    RealScalarFunction zero=RealScalarFunction::constant(2,0.0);
    RealScalarFunction one=RealScalarFunction::constant(2,1.0);
    RealScalarFunction two=RealScalarFunction::constant(2,2.0);
    RealScalarFunction x=RealScalarFunction::coordinate(2,0);
    RealScalarFunction y=RealScalarFunction::coordinate(2,1);

    MonolithicHybridAutomaton ha;
    ha.new_mode(location1,(one,one));
    ha.new_mode(location2,(one,-one));
    ha.new_mode(location3,(one,two));
    ha.new_mode(location4,(one,zero));
    ha.new_transition(location1,event1,location2,(x,y-1),(x-0.53125),urgent);
    ha.new_transition(location2,event2,location3,(x,y+1),(x-1.53125),urgent);

    // Define a bounding box for the evolution
    std::cout<<"making bounding_box"<<std::endl;
    Box bounding_box=Box(2, -1.0,3.0, -2.0,2.0);
    std::cout<<"bounding_box="<<bounding_box<<"\n"<<std::endl;
    //Box eps_bounding_box=bounding_box.neighbourhood(0.1);

    // Define the initial cell
    //Box initial_box=make_box("[-0.0625,+0.0625]x[-0.0625,+0.0625]");
    Box initial_box=make_box("[0.001,0.002]x[0.001,0.002]");
    ARIADNE_TEST_PRINT(initial_box);
    GridTreeSet approx_tree_set=outer_approximation(initial_box,grid,depth);
    GridCell initial_cell=*approx_tree_set.begin();
    HybridGridCell hybrid_initial_cell(location1,initial_cell);
    ARIADNE_TEST_PRINT(hybrid_initial_cell);
    HybridBox hybrid_initial_set=hybrid_initial_cell.box();
    ARIADNE_TEST_PRINT(hybrid_initial_set);
    //[1.00098:1.00122],
    HybridTime htime(time,steps);
    ARIADNE_TEST_PRINT(htime);

    // Compute the reachable sets
    cout << "Computing evolution... " << flush;
    // evolver.verbosity=1;
    Orbit<EnclosureType> evolve_orbit
        = evolver.orbit(ha,hybrid_initial_set,htime,UPPER_SEMANTICS);
    cout << "done." << endl;

    ARIADNE_TEST_PRINT(evolve_orbit);

    cout << "Extracting grid... " << flush;
    HybridGrid hagrid(ha.state_space());
    cout << "done." << endl;

    // Compute the reachable sets
    cout << "Computing discretised evolution... " << flush;
    Orbit<HybridGridCell> discrete_orbit
        = discrete_evolver.evolution(ha,hybrid_initial_cell,htime,depth,UPPER_SEMANTICS);
    cout << "done." << endl;

    ContinuousEnclosureType const& initial_set=evolve_orbit.initial().continuous_state_set();
    ListSet<ContinuousEnclosureType> const& reach_set1=evolve_orbit.reach()[location1];
    ListSet<ContinuousEnclosureType> const& reach_set2=evolve_orbit.reach()[location2];
    ListSet<ContinuousEnclosureType> const& reach_set3=evolve_orbit.reach()[location3];
    ListSet<ContinuousEnclosureType> const& final_set=evolve_orbit.final()[location3];

    GridTreeSet const& reach_cells1=discrete_orbit.reach()[location1];
    GridTreeSet const& reach_cells2=discrete_orbit.reach()[location2];
    GridTreeSet const& reach_cells3=discrete_orbit.reach()[location3];
    GridTreeSet const& final_cells=discrete_orbit.final()[location3];


    ARIADNE_TEST_PRINT(initial_set);
    ARIADNE_TEST_PRINT(initial_cell);
    ARIADNE_TEST_PRINT(reach_set1);
    ARIADNE_TEST_PRINT(reach_cells1);
    ARIADNE_TEST_PRINT(final_set);
    ARIADNE_TEST_PRINT(final_cells);


    cout << "Plotting... " << flush;

    // Plot the intial, evolve and reach sets
    {
        Figure fig;
        fig.set_bounding_box(bounding_box);
        fig << line_style(true);
        fig << line_style(true);
        fig << fill_colour(cyan) << reach_cells1;
        fig << fill_colour(cyan) << reach_cells2;
        fig << fill_colour(cyan) << reach_cells3;
        fig << fill_colour(blue) << initial_cell;
        fig << fill_colour(magenta) << final_cells;
        fig << fill_colour(cyan) << reach_set1;
        fig << fill_colour(cyan) << reach_set2;
        fig << fill_colour(cyan) << reach_set3;
        fig << fill_colour(blue) << initial_set;
        fig << fill_colour(magenta) << final_set;
        fig.write("test_hybrid_discretiser");
    }

    cout << "done." << endl;

    ARIADNE_TEST_PRINT(evolve_orbit.reach()[location4]);
    ARIADNE_TEST_PRINT(discrete_orbit.reach()[location4]);
}

