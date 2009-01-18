/***************************************************************************
 *            test_continuous_evolution.cc
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
#include "function.h"
#include "taylor_variable.h"
#include "taylor_set.h"
#include "taylor_function.h"
#include "box.h"
#include "zonotope.h"
#include "list_set.h"
#include "evolution_parameters.h"
#include "taylor_calculus.h"
#include "vector_field_evolver.h"
#include "orbit.h"
#include "graphics.h"
#include "logging.h"

#include "models.h"

#include "test.h"

using namespace Ariadne;
using namespace std;
using Models::Henon;

class TestContinuousEvolution
{
  public:
    void test() const;
    void simple_test() const;
};

int main() 
{
    TestContinuousEvolution().simple_test();
    TestContinuousEvolution().test();
    return ARIADNE_TEST_FAILURES;
}

void TestContinuousEvolution::simple_test() const
{
    std::cout <<std::endl; 

    TaylorCalculus calculus;
    AffineFunction vector_field(Matrix<Float>(1,1,0.0),Vector<Float>(1,1.0));
    Vector<Interval> initial_box(1,Interval(-0.01,0.01));
    TaylorSet initial_set(initial_box);
    double step_size;
    Box bounding_box(2);
    make_lpair(step_size,bounding_box)=calculus.flow_bounds(vector_field,initial_box,0.5,1.0);
    step_size=0.25;
    bounding_box=Vector<Interval>(1,Interval(-0.5,1.0));
    TaylorFunction vector_field_expansion=calculus.map_model(vector_field,Vector<Interval>(1,Interval(-1,+1)));
    TaylorFunction vector_field_model=calculus.map_model(vector_field,bounding_box);
    TaylorFunction flow_model=calculus.flow_model(vector_field,initial_box,step_size,bounding_box);
    TaylorSet evolve_set=calculus.integration_step(flow_model,initial_set,step_size);
    TaylorSet reach_set=calculus.reachability_step(flow_model,initial_set,0.0,step_size);
    std::cout <<"vector_field="<<vector_field<<std::endl;
    std::cout <<"vector_field_expansion="<<vector_field_expansion<<std::endl;
    std::cout <<"step_size="<<step_size<<std::endl;
    std::cout <<"bounding_box="<<bounding_box<<std::endl;
    std::cout <<"vector_field_model="<<vector_field_model<<std::endl;
    std::cout <<"vector_field_model_expansion="<<expansion(vector_field_model.variables(),vector_field_model.domain())<<std::endl;
    std::cout <<"flow_model="<<flow_model<<std::endl;
    std::cout <<"flow_model_expansion="<<expansion(flow_model.variables(),flow_model.domain())<<std::endl;
    std::cout <<"evolve_set="<<evolve_set<<std::endl;
    std::cout <<"reach_set="<<reach_set<<std::endl;
    //std::cout <<"="<<<<std::endl;
}

    

void TestContinuousEvolution::test() const
{
    // cout << __PRETTY_FUNCTION__ << endl;

    typedef TaylorSet EnclosureType;

    // Set up the evolution parameters and grid
    Float time(5.0);
    Float step_size(0.125);
    Float enclosure_radius(0.25);
    
    EvolutionParameters parameters;
    parameters.maximum_enclosure_radius=enclosure_radius;
    parameters.maximum_step_size=step_size;

    // Set up the evaluators
    VectorFieldEvolver evolver(parameters);
  
    // Define the initial box
    Box initial_box(2); 
    initial_box[0]=Interval(1.01,1.02);
    initial_box[1]=Interval(0.51,0.52);

    // cout << "initial_box=" << initial_box << endl;

    // Set up the vector field
    Float mu=0.5;
    Vector<Float> p(1); p[0]=mu;
    Function<VanDerPol> vdp(p);
    // cout << "van_der_pol_function=" << vdp << endl;
    // cout << "van_der_pol_function.parameters()=" << vdp.parameters() << endl;

    //Function evaluation sanity check
    // cout << "vdp.evaluate(" << initial_box << ") " << flush; // cout << " = " << vdp.evaluate(initial_box) << endl;
    // cout << "vdp.jacobian(" << initial_box << ") = " << vdp.jacobian(initial_box) << endl;
    // cout << endl;
  
    AffineFunction aff(Matrix<Float>(2,2,0.,0.,1.,0.),Vector<Float>(2,1.,0.0));

    // Make a hybrid automaton for the Van der Pol equation
    //VectorField vanderpol(aff);
    VectorField vanderpol(vdp);


    // Over-approximate the initial set by a grid cell
    EnclosureType initial_set(initial_box);
    // cout << "initial_set=" << initial_set << endl << endl;
  
    Semantics semantics=UPPER_SEMANTICS;

    // Compute the reachable sets
    ListSet<EnclosureType> evolve_set,reach_set;
    Orbit<EnclosureType> orbit = evolver.orbit(vanderpol,initial_set,time,semantics);

    cout << "\norbit=\n" << orbit << endl << endl;

    // Print the intial, evolve and reach sets
    // cout << "Plotting sets" << endl;
    // cout << "evolve_set=" << hybrid_evolve_set << endl;
    // cout << "reach_set=" << hybrid_reach_set << endl;
    Figure fig;
    fig << line_style(true) << fill_colour(cyan) << orbit.reach();
    fig << fill_colour(magenta) << orbit.intermediate();
    fig << fill_colour(red) << orbit.final();
    fig << fill_colour(blue) << initial_set;
    fig.write("test_continuous_evolution-vdp");

}
