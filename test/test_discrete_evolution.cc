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
#include "sparse_differential.h"
#include "differential_vector.h"
#include "function.h"
#include "approximate_taylor_model.h"
#include "box.h"
#include "zonotope.h"
#include "list_set.h"
#include "evolution_parameters.h"
#include "hybrid_evolver.h"
#include "graphics.h"
#include "logging.h"

#include "models.h"

#include "test.h"

using namespace Ariadne;
using namespace std;
using Models::Henon;

class TestDiscreteEvolver
{
  public:
    void test() const;
};

int main() 
{
    TestDiscreteEvolver().test();
    return ARIADNE_TEST_FAILURES;
}

void TestDiscreteEvolver::test() const
{
    cout << __PRETTY_FUNCTION__ << endl;

    typedef ApproximateTaylorModel EnclosureType;
    typedef std::pair<DiscreteState,ApproximateTaylorModel> HybridEnclosureType;

    // Set up the evolution parameters and grid
    Float time(6.0);
    uint steps(5);
    Float step_size(0.0625);
    Float enclosure_radius(0.25);
    
    EvolutionParameters parameters;
    parameters.maximum_enclosure_radius=enclosure_radius;
    parameters.maximum_step_size=step_size;

    // Set up the evaluators
    HybridEvolver evolver(parameters);
  
    // Define the initial box
    Box initial_box(2); 
    initial_box[0]=Interval(1.01,1.02);
    initial_box[1]=Interval(0.51,0.52);

    cout << "initial_box=" << initial_box << endl;

    // Set up the vector field
    Vector<Float> p(2); p[0]=1.5; p[1]=0.375;
    Function<Henon> h(p);
    cout << "henon_function=" << h << endl;
    cout << "henon_function.parameters()=" << h.parameters() << endl;

    //Function evaluation sanity check
    Vector<Float> x(2); x[0]=0.5; x[1]=0.25; 
    Vector<Float> hx(2); hx[0]=p[0]-x[0]*x[0]+x[1]*p[1]; hx[1]=x[0];
    ARIADNE_TEST_EQUAL(h.evaluate(x),hx);
    Matrix<Float> dhx(2,2); dhx[0][0]=-2*x[0]; dhx[0][1]=p[1]; dhx[1][0]=1.0;
    ARIADNE_TEST_EQUAL(h.jacobian(x),dhx);
 

    //Function evaluation sanity check
    cout << "h.evaluate(" << initial_box << ") " << flush; cout << " = " << h.evaluate(initial_box) << endl;
    cout << "h.jacobian(" << initial_box << ") = " << h.jacobian(initial_box) << endl;

  
    // Make a hybrid automaton for the Henon function
    HybridAutomaton henon("Henon");
    DiscreteState location(42);
    henon.new_mode(location,ConstantFunction(Vector<Float>(2),2));


    // Over-approximate the initial set by a grid cell
    EnclosureType initial_set(initial_box,IdentityFunction(2),4,1);
    cout << "initial_set=" << initial_set << endl << endl;
    HybridEnclosureType initial_hybrid_set(location,initial_set);
    HybridTime hybrid_time(time,steps);

  
    // Compute the reachable sets
    ListSet<HybridEnclosureType> hybrid_evolve_set,hybrid_reach_set;
    hybrid_evolve_set = evolver.evolve(henon,initial_hybrid_set,hybrid_time);
    //cout << "evolve_set=" << hybrid_evolve_set << endl;
    hybrid_reach_set = evolver.reach(henon,initial_hybrid_set,hybrid_time);
    //cout << "reach_set=" << hybrid_reach_set << endl;
  
    // Print the intial, evolve and reach sets
    Figure fig;
    fig << line_style(true) << fill_colour(cyan) << hybrid_reach_set;
    fig << fill_colour(yellow) << hybrid_evolve_set;
    fig << fill_colour(blue) << initial_set;
    fig.write("test_discrete_evolution-henon");
}
