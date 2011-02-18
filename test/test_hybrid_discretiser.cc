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

#include "tuple.h"
#include "vector.h"
#include "matrix.h"
#include "function.h"
#include "taylor_set.h"
#include "zonotope.h"
#include "list_set.h"
#include "grid_set.h"
#include "graphics.h"
#include "hybrid_time.h"
#include "hybrid_set.h"
#include "hybrid_automaton.h"
#include "hybrid_evolver.h"
#include "hybrid_orbit.h"
#include "hybrid_discretiser.h"

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
    std::cerr<<"SKIPPED "; return 1u;
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

    cout << __PRETTY_FUNCTION__ << endl;

    // Set up the evolution parameters and grid
    Real time(1.0);
    uint steps(6);
    double maximum_step_size(0.125);
    int depth=8;
    DiscreteLocation location(1);
    DiscreteEvent event(1);

    EvolutionParameters parameters;
    parameters.maximum_step_size=maximum_step_size;
    Grid grid(2);

    // Set up the evaluators
    EvolverType evolver(parameters);
    HybridDiscretiser< EnclosureType > discrete_evolver(evolver);


    // Set up the vector field
    Real a=1.5; Real b=0.375;
    RealScalarFunction zero=RealScalarFunction::constant(2,0.0);
    RealScalarFunction one=RealScalarFunction::constant(2,1.0);
    RealScalarFunction x=RealScalarFunction::coordinate(2,0);
    RealScalarFunction y=RealScalarFunction::coordinate(2,1);

    MonolithicHybridAutomaton ha("Decay");
    ha.new_mode(location,(one,-y));
    ha.new_transition(location,event,location,(x-1,y),(x-1),urgent);

    // Define a bounding box for the evolution
    std::cout<<"making bounding_box"<<std::endl;
    Box bounding_box=make_box("[-4,4]x[-4,4]") ;
    std::cout<<"bounding_box="<<bounding_box<<"\n"<<std::endl;
    //Box eps_bounding_box=bounding_box.neighbourhood(0.1);

    // Define the initial cell
    Box initial_box=make_box("[1.0001,1.0002]x[0.5001,0.5002]");
    ARIADNE_TEST_PRINT(initial_box);
    GridTreeSet approx_tree_set=outer_approximation(initial_box,grid,depth);
    GridCell initial_cell=*approx_tree_set.begin();
    HybridGridCell hybrid_initial_cell(location,initial_cell);
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
        = evolver.orbit(ha,EnclosureType(hybrid_initial_set),htime,UPPER_SEMANTICS);
    cout << "done." << endl;

    ARIADNE_TEST_PRINT(evolve_orbit);

    cout << "Extracting grid... " << flush;
    HybridGrid hagrid=ha.grid();
    cout << "done." << endl;

    // Compute the reachable sets
    cout << "Computing discretised evolution... " << flush;
    Orbit<HybridGridCell> discrete_orbit
        = discrete_evolver.evolution(ha,hybrid_initial_cell,htime,depth,UPPER_SEMANTICS);
    cout << "done." << endl;

    ContinuousEnclosureType const& initial_set=evolve_orbit.initial().continuous_state_set();
    ListSet<ContinuousEnclosureType> const& reach_set=evolve_orbit.reach()[location];
    ListSet<ContinuousEnclosureType> const& intermediate_set=evolve_orbit.intermediate()[location];
    ListSet<ContinuousEnclosureType> const& final_set=evolve_orbit.final()[location];

    GridTreeSet const& reach_cells=discrete_orbit.reach()[location];
    GridTreeSet const& intermediate_cells=discrete_orbit.intermediate()[location];
    GridTreeSet const& final_cells=discrete_orbit.final()[location];


    ARIADNE_TEST_PRINT(initial_set);
    ARIADNE_TEST_PRINT(initial_cell);
    ARIADNE_TEST_PRINT(reach_set);
    ARIADNE_TEST_PRINT(reach_cells);
    ARIADNE_TEST_PRINT(intermediate_set);
    ARIADNE_TEST_PRINT(intermediate_cells);
    ARIADNE_TEST_PRINT(final_set);
    ARIADNE_TEST_PRINT(final_cells);



    // Plot the intial, evolve and reach sets
    {
        Figure fig;
        fig.set_bounding_box(Box(2,Interval(-3,3)));
        fig << line_style(true);
        fig << line_style(true);
        fig << fill_colour(cyan) << reach_cells;
        fig << fill_colour(magenta) << intermediate_cells;
        fig << fill_colour(yellow) << initial_cell;
        fig << fill_colour(green) << final_cells;
        fig.write("test_discrete_evolver-hybrid-cells");
    }

    // Plot the intial, evolve and reach sets
    {
        Figure fig;
        fig.set_bounding_box(Box(2,Interval(-3,3)));
        fig << line_style(true);
        fig << fill_colour(cyan) << reach_set;
        fig << fill_colour(magenta) << intermediate_set;
        fig << fill_colour(yellow) << initial_set;
        fig << fill_colour(green) << final_set;
        fig.write("test_discrete_evolver-hybrid-sets");
    }

}

