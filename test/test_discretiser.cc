/***************************************************************************
 *            test_discretiser.cc
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
#include "integrator.h"
#include "map_evolver.h"
#include "vector_field_evolver.h"
#include "orbit.h"
#include "discretiser.h"
#include "graphics.h"

#include "logging.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

class TestDiscretiser
{
  public:
    void test() const;
  private:
    void test_discrete_time() const;
    void test_continuous_time() const;
};

int main()
{
    std::cerr<<"SKIPPED "; return 1u;
    TestDiscretiser().test();
    return ARIADNE_TEST_FAILURES;
}

void TestDiscretiser::test() const
{
    ARIADNE_TEST_CALL(test_discrete_time());
    ARIADNE_TEST_CALL(test_continuous_time());
}


void TestDiscretiser::test_discrete_time() const
{
    typedef MapEvolver::EnclosureType EnclosureType;

    cout << __PRETTY_FUNCTION__ << endl;

    // Set up the evolution parameters and grid
    uint steps(6);
    double maximum_step_size(0.125);
    int depth=8;

    EvolutionParameters parameters;
    parameters.maximum_step_size=maximum_step_size;
    Grid grid(2);

    // Set up the evaluators
    MapEvolver evolver(parameters);

    Discretiser< IteratedMap, EnclosureType > discrete_evolver(evolver);


    // Set up the vector field
    Real a=1.5; Real b=0.375;
    RealScalarFunction x=RealScalarFunction::coordinate(2,0);
    RealScalarFunction y=RealScalarFunction::coordinate(2,1);
    RealVectorFunction henon=join(a-x*x+b*y,x);
    cout << "henon=" << henon << endl;
    IteratedMap system(henon);

    // Define a bounding box for the evolution
    std::cout<<"making bounding_box"<<std::endl;
    Box bounding_box=make_box("[-4,4]x[-4,4]") ;
    std::cout<<"bounding_box="<<bounding_box<<"\n"<<std::endl;

    // Define the initial cell
    Box box=make_box("[1.001,1.002]x[0.501,0.502]");
    GridTreeSet approx_tree_set=outer_approximation(box,grid,depth);
    GridCell initial_cell=*approx_tree_set.begin();
    cout << "initial_cell=" << initial_cell << endl << endl;
    Box initial_box=initial_cell.box();
    cout << "initial_box=" << initial_box << endl << endl;
    //[1.00098:1.00122],
    cout << "steps=" << steps << endl << endl;

    // Compute the reachable sets
    cout << "Computing evolution... " << flush;
    Orbit<EnclosureType> evolve_orbit
        = evolver.orbit(system,EnclosureType(initial_box,Sweeper()),steps,UPPER_SEMANTICS);
    cout << "done." << endl;

    EnclosureType const& initial_set=evolve_orbit.initial();
    ListSet<EnclosureType> const& reach_set=evolve_orbit.reach();
    ListSet<EnclosureType> const& intermediate_set=evolve_orbit.intermediate();
    ListSet<EnclosureType> const& final_set=evolve_orbit.final();

    // Compute the reachable sets
    cout << "Computing discretised evolution... " << flush;
    Orbit<GridCell> discrete_orbit
        = discrete_evolver.upper_evolution(system,initial_cell,steps,grid,depth);
    cout << "done." << endl;

    GridTreeSet const& reach_cells=discrete_orbit.reach();
    GridTreeSet const& intermediate_cells=discrete_orbit.intermediate();
    GridTreeSet const& final_cells=discrete_orbit.final();

    cout << "initial_set=" << initial_set.bounding_box() << endl << endl;
    cout << "initial_cell=" << initial_cell.box() << endl << endl;
    cout << "reach_set=" << reach_set << endl << endl;
    cout << "reach_cells=" << reach_cells << endl << endl;
    cout << "intermediate_set=" << intermediate_set << endl << endl;
    cout << "intermediate_cells=" << intermediate_cells << endl << endl;
    cout << "final_set=" << final_set << endl << endl;
    cout << "final_cells=" << final_cells << endl << endl;

    // Plot the intial, evolve and reach sets
    {
        Figure fig;
        fig.set_bounding_box(Box(2,Interval(-3,3)));
        fig << line_style(true);
        fig << line_style(true);
        fig << fill_colour(cyan) << reach_cells;
        fig << fill_colour(yellow) << initial_cell;
        fig << fill_colour(green) << final_cells;
        fig.write("test_discretiser-henon-cells");
    }

    // Plot the intial, evolve and reach sets
    {
        Figure fig;
        fig.set_bounding_box(Box(2,Interval(-3,3)));
        fig << line_style(true);
        fig << fill_colour(cyan) << reach_set;
        fig << fill_colour(yellow) << initial_set;
        fig << fill_colour(green) << final_set;
        fig.write("test_discretiser-henon-sets");
    }
}

void TestDiscretiser::test_continuous_time() const
{
    typedef TaylorConstrainedImageSet EnclosureType;

    cout << __PRETTY_FUNCTION__ << endl;

    // Set up the evolution parameters and grid
    Real time(1.0);
    double maximum_step_size(0.125);
    int depth=8;

    EvolutionParameters parameters;
    parameters.maximum_step_size=maximum_step_size;
    Grid grid(2);

    // Set up the evaluators
    TaylorIntegrator integrator(4,1e-4);
    VectorFieldEvolver evolver(parameters,integrator);
    Discretiser< VectorField, EnclosureType > discretiser(evolver);


    // Set up the vector field
    Real mu=0.865; Vector<Real> p(1); p[0]=mu;
    RealVectorFunction x=RealVectorFunction::identity(2);
    RealVectorFunction vdp = ( x[1]+Real(0), p[0]*(1-x[0]*x[0])*x[1]-x[0] );
    cout << "vdp=" << vdp << endl;
    VectorField system(vdp);

    // Define a bounding box for the evolution
    Box bounding_box=make_box("[-4,4]x[-4,4]") ;
    //Box eps_bounding_box=bounding_box.neighbourhood(0.1);

    // Define the initial cell
    Box box=make_box("[1.01,1.02]x[0.51,0.52]");
    cout << "box=" << box << endl;
    GridTreeSet approx_tree_set=outer_approximation(box,grid,depth);
    GridCell initial_cell=*approx_tree_set.begin();
    Box initial_box=initial_cell.box();
    cout << "initial_box=" << initial_box << endl << endl;

    // Compute the reachable sets
    cout << "Computing evolution... " << flush;
    Orbit<EnclosureType> evolve_orbit
        = evolver.orbit(system,EnclosureType(initial_box,Sweeper()),time,UPPER_SEMANTICS);
    cout << "done." << endl;
    EnclosureType const& initial_set=evolve_orbit.initial();
    ListSet<EnclosureType> const& reach_set=evolve_orbit.reach();
    ListSet<EnclosureType> const& intermediate_set=evolve_orbit.intermediate();
    ListSet<EnclosureType> const& final_set=evolve_orbit.final();

    // Compute the reachable sets
    cout << "Computing discretised evolution... " << flush;
    Orbit<GridCell> discrete_orbit
        = discretiser.upper_evolution(system,initial_cell,time,grid,depth);
    cout << "done." << endl;
    GridTreeSet const& reach_cells=discrete_orbit.reach();
    GridTreeSet const& intermediate_cells=discrete_orbit.intermediate();
    GridTreeSet const& final_cells=discrete_orbit.final();

    cout << "initial_set_bounding_box=" << initial_set.bounding_box() << endl << endl;
    cout << "initial_cell=" << initial_cell.box() << endl << endl;
    cout << "final_set=" << final_set << endl << endl;
    cout << "final_cells=" << final_cells << endl << endl;

    // Plot the intial, evolve and reach sets
    {
        Figure fig;
        fig.set_bounding_box(Box(2,Interval(-3,3)));
        fig << line_style(true);
        fig << line_style(true);
        fig << fill_colour(cyan) << reach_cells;
        fig << fill_colour(magenta) << intermediate_cells;
        fig << fill_colour(blue) << initial_cell;
        fig << fill_colour(blue) << final_cells;
        fig.write("test_discretiser-vdp-cells");
    }

    // Plot the intial, evolve and reach sets
    {
        Figure fig;
        fig.set_bounding_box(Box(2,Interval(-3,3)));
        fig << line_style(true);
        fig << fill_colour(cyan) << reach_set;
        fig << fill_colour(magenta) << intermediate_set;
        fig << fill_colour(blue) << initial_set;
        fig << fill_colour(blue) << final_set;
        fig.write("test_discretiser-vdp-sets");
    }


}


