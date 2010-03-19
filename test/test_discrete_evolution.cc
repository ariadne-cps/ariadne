/***************************************************************************
 *            test_discrete_evolution.cc
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

#include <fstream>
#include <iostream>

#include "tuple.h"
#include "vector.h"
#include "matrix.h"
#include "taylor_model.h"
#include "differential.h"
#include "function.h"
#include "taylor_function.h"
#include "taylor_set.h"
#include "box.h"
#include "zonotope.h"
#include "list_set.h"
#include "evolution_parameters.h"
#include "map.h"
#include "map_evolver.h"
#include "graphics.h"
#include "logging.h"

#include "models.h"

#include "test.h"

using namespace Ariadne;
using namespace std;
using Models::Henon;

class TestMapEvolver
{
  public:
    void test() const;
};

int main()
{
    TestMapEvolver().test();
    return ARIADNE_TEST_FAILURES;
}

void TestMapEvolver::test() const
{
    cout << __PRETTY_FUNCTION__ << endl;

    typedef TaylorImageSet EnclosureType;

    // Set up the evolution parameters and grid
    IteratedMap::TimeType time(3);
    Float enclosure_radius(0.25);

    EvolutionParameters parameters;
    parameters.maximum_enclosure_radius=enclosure_radius;

    // Set up the evaluators
    MapEvolver evolver(parameters);

    // Define the initial box
    Box initial_box(2);
    initial_box[0]=Interval(1.01,1.03);
    initial_box[1]=Interval(0.51,0.53);

    cout << "initial_box=" << initial_box << endl;

    // Set up the vector field
    Vector<Float> p(2); p[0]=1.5; p[1]=0.375;
    VectorUserFunction<Henon> henon(p);
    cout << "henon_function=" << henon << endl;
    cout << "henon_function.parameters()=" << henon.parameters() << endl;

    //VectorUserFunction evaluation sanity check
    Vector<Float> x(2); x[0]=0.5; x[1]=0.25;
    Vector<Float> hx(2); hx[0]=p[0]-x[0]*x[0]+x[1]*p[1]; hx[1]=x[0];
    ARIADNE_TEST_EQUAL(henon.evaluate(x),hx);
    Matrix<Float> dhx(2,2); dhx[0][0]=-2*x[0]; dhx[0][1]=p[1]; dhx[1][0]=1.0;
    ARIADNE_TEST_EQUAL(henon.jacobian(x),dhx);


    //VectorUserFunction evaluation sanity check
    cout << "henon.evaluate(" << initial_box << ") " << flush; cout << " = " << henon.evaluate(initial_box) << endl;
    cout << "henon.jacobian(" << initial_box << ") = " << henon.jacobian(initial_box) << endl;



    // Over-approximate the initial set by a grid cell
    EnclosureType initial_set(initial_box);
    cout << "initial_set=" << initial_set << endl << endl << endl;


    // Compute the reachable sets
    ListSet<EnclosureType> evolve_set,reach_set;
    //evolve_set = evolver.evolve(henon,initial_set,time);
    reach_set = evolver.reach(henon,initial_set,time);
    cout << "initial_bounding_box=" << initial_set.bounding_box() << endl;
    //cout << "evolve_bounding_boxes=" << evolve_set.bounding_boxes() << endl;
    cout << "reach_bounding_boxes=" << reach_set.bounding_boxes() << endl;

    // Print the intial, evolve and reach sets
    Figure fig;
    fig << line_style(true) << fill_colour(cyan) << reach_set;
    fig << fill_colour(yellow) << evolve_set;
    fig << fill_colour(blue) << initial_set;
    fig.write("test_discrete_evolution-henon");
}
