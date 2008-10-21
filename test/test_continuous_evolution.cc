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

class TestContinuousEvolution
{
 public:
  void test() const;
};

int main() 
{
  TestContinuousEvolution().test();
  return ARIADNE_TEST_FAILURES;
}

void TestContinuousEvolution::test() const
{
  cout << __PRETTY_FUNCTION__ << endl;

  typedef ApproximateTaylorModel EnclosureType;
  typedef std::pair<DiscreteState,ApproximateTaylorModel> HybridEnclosureType;

  // Set up the evolution parameters and grid
  Float time(4.0);
  Float step_size(0.0625);
  Float enclosure_radius(0.25);
    
  EvolutionParameters parameters;
  parameters.maximum_enclosure_radius=enclosure_radius;
  parameters.maximum_step_size=step_size;

  // Set up the evaluators
  HybridEvolver evolver;
  
  // Define the initial box
  Box initial_box(2); 
  initial_box[0]=Interval(1.01,1.02);
  initial_box[1]=Interval(0.51,0.52);

  cout << "initial_box=" << initial_box << endl;

  // Set up the vector field
  Float mu=0.5;
  Vector<Float> p(1); p[0]=mu;
  Function<VanDerPol> vdp(p);
  cout << "van_der_pol_function=" << vdp << endl;
  cout << "van_der_pol_function.parameters()=" << vdp.parameters() << endl;

  //Function evaluation sanity check
  cout << "vdp.evaluate(" << initial_box << ") " << flush; cout << " = " << vdp.evaluate(initial_box) << endl;
  cout << "vdp.jacobian(" << initial_box << ") = " << vdp.jacobian(initial_box) << endl;
  cout << endl;
  
  // Make a hybrid automaton for the Henon function
  DiscreteState location(23);

  // Make a hybrid automaton for the Henon function
  HybridAutomaton vanderpol("Van der Pol");
  vanderpol.new_mode(location,vdp);


  // Over-approximate the initial set by a grid cell
  Vector<Interval> unit_box(2,Interval(-1,1));
  EnclosureType initial_set=ApproximateTaylorModel(unit_box,ScalingFunction(initial_box),4,1);
  cout << "initial_set=" << initial_set << endl << endl;
  HybridEnclosureType initial_hybrid_set(location,initial_set);
  HybridTime hybrid_time(1.0,5);

  
  // Compute the reachable sets
  ListSet<HybridEnclosureType> hybrid_evolve_set,hybrid_reach_set;
  hybrid_evolve_set = evolver.evolve(vanderpol,initial_hybrid_set,hybrid_time);
  cout << "evolve_set=" << hybrid_evolve_set << endl;
  hybrid_reach_set = evolver.reach(vanderpol,initial_hybrid_set,hybrid_time);
  cout << "reach_set=" << hybrid_reach_set << endl;
  
  //cout << "evolve_set=" << hybrid_evolve_set[location] << endl;
  //cout << "reach_set=" << hybrid_reach_set[location] << endl;

  // Print the intial, evolve and reach sets
  cout << "Plotting sets" << endl;
  std::ofstream ofs("test_continuous_evolution-vdp");
  //Graphic plot(ofstream);
  Graphic fig;
  //fig.set_bounding_box(Box(hybrid_reach_set[0].second.range()));
  //fig << line_style(true) << fill_colour(cyan) << hybrid_reach_set[0].second;
  //fig << fill_colour(yellow) << hybrid_evolve_set[0].second;
  //fig << fill_colour(blue) << initial_set;
  Box bx1,bx2;
  bx1=hybrid_reach_set[0].second.range();
  bx2=hybrid_evolve_set[0].second.range();
  fig.plot(bx1);
  fig.plot(bx2);
  //fig << bx2;
  fig.display();
  ofs.close();

}
