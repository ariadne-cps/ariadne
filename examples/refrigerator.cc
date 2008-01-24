/***************************************************************************
 * refrig_sync.cc
 * What: Case study of refrigiration syncronization.
 * Who:  Jesper A. Larsen (Aalborg University/DK & CWI/NL) jal@es.aau.dk
 * When: 20070928
 *
 ****************************************************************************/


#include <iostream>
#include <fstream>
#include <string>

#include "ariadne.h"
#include "numeric/float64.h"
#include "geometry/linear_constraint.h"
#include "geometry/set_interface.h"
#include "geometry/set_reference.h"
#include "geometry/hybrid_set.h"
#include "geometry/box.h"
#include "geometry/polyhedron.h"
#include "geometry/polyhedral_set.h"
#include "geometry/constraint_set.h"
#include "function/function_interface.h"
#include "system/affine_map.h"
#include "system/affine_vector_field.h"
#include "system/hybrid_automaton.h"
//#include "system/constraint_based_hybrid_automaton.h"
#include "evaluation/evolution_parameters.h"
#include "evaluation/standard_applicator.h"
#include "evaluation/lohner_integrator.h"
#include "evaluation/affine_integrator.h"
#include "evaluation/set_based_hybrid_evolver.h"
//#include "evaluation/constraint_based_hybrid_evolver.h"
#include "output/epsstream.h"
#include "output/logging.h"
#define alpha1 0.15
#define alpha2 0.15

using namespace Ariadne;
using namespace Ariadne::Function;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;
using namespace Ariadne::Output;
using namespace Ariadne::Combinatoric;
using namespace std;

template<class R> int refrig();
  
int main() {
  return refrig<Float64>();
}

template<class R>
int refrig() 
{
  set_hybrid_evolver_verbosity(4);
  typedef Zonotope<R> BS;
  
  Box<R> domain("[-1,1]x[-1,1]"); //Setup domain
  cout << "domain=" << domain << endl;

  AffineVectorField<R> dynamic1(Matrix<R>("[0,0;0,0]"),Vector<R>("[1,1]")); //First dyn
  cout << "1st dynamic=" << dynamic1 << endl;

  double A[4]={-2*alpha1,0,-2*alpha2,0};
  double v[2]={-1,1};
  AffineVectorField<R> dynamic2(Matrix<R>(2,2,A),Vector<R>(2,v));
  cout << "2nd dynamic=" << dynamic2 << endl;

	double B[4]={0,-2*alpha1,0,-2*alpha2};
  AffineVectorField<R> dynamic3(Matrix<R>(2,2,B),Vector<R>("[1,-1]")); //Third dyn

	double C[4]={-2*alpha1,-2*alpha2,-2*alpha1,-2*alpha2};
  AffineVectorField<R> dynamic4(Matrix<R>(2,2,C),Vector<R>("[-1,-1]")); //Fourth dyn

  Vector<R> e1=Vector<R>::unit(2,1);
  LinearConstraint<R> constraint1(e1,1);
  RectangularSet<R> inv(Box<R>("[-1,1]x[-1,1]")); // All sets have the same invariants.

  cout << "invariant=" << inv << endl;

  RectangularSet<R> guard(Box<R>("[0,0]x[1,1]"));

  AffineMap<R> reset(Matrix<R>("[1,0;0,1]"),Vector<R>("[0,0]"));
  cout << "reset=" << reset << endl;
  
  HybridAutomaton<R> automaton("Refridguration desyncronization");
  // Building continues modes
  DiscreteState state1(1);
  DiscreteState state2(2);
  DiscreteState state3(3);
  DiscreteState state4(4);
  cout << "states={" << state1 << " " << state2 << " " << state3 << " " << state4 << "}" << endl;

  //  DiscretState state1=mode1_id;
  const DiscreteMode<R>& mode1=automaton.new_mode(state1,dynamic1,inv);
  const DiscreteMode<R>& mode2=automaton.new_mode(state2,dynamic2,inv);
  const DiscreteMode<R>& mode3=automaton.new_mode(state3,dynamic3,inv);
  const DiscreteMode<R>& mode4=automaton.new_mode(state4,dynamic4,inv);
  cout << "modes={" << mode1 << " " << mode2 << " " << mode3 << " " << mode4 << "}" << endl;

  // Form the forced transitions 
  DiscreteEvent event1(5);
  DiscreteEvent event2(6);
  DiscreteEvent event3(7);
  DiscreteEvent event4(8);
  DiscreteEvent event5(9);
  DiscreteEvent event6(10);
  cout << "events={" << event1 << " " << event2 << " " << event3 << " " << event4 << " " << event5 << " " << event6 << "}" << endl;

  //const ConstraintBasedDiscreteTransition<R>& transition1=automaton.new_forced_transition(event1,state1,state2,reset,guard);

  PolyhedralSet<R> invariant1(constraint1.polyhedron());
  cout << "invariant1=" << invariant1 << endl;
  const DiscreteTransition<R>& transition1=automaton.new_transition(event1,state1,state2,reset,invariant1);
  cout << "transition1=" << transition1 << endl;

  EvolutionParameters<R> parameters;
  //cout << "parameters.minimum_step_size()=" << flush; cout << Rational(parameters.minimum_step_size()) << endl;
  parameters.set_maximum_step_size(0.125);
  parameters.set_lock_to_grid_time(0.5);
  parameters.set_maximum_basic_set_radius(0.25);
  parameters.set_grid_length(0.125);
  cout << "parameters=" << parameters << endl;


  SetBasedHybridEvolver<BS> hybrid_evolver(parameters);
  
  RectangularSet<R> bounding_box(domain.neighbourhood(1));
  cout << "bounding_box=" << bounding_box << endl;
  
  RectangularSet<R> initial_box(Box<R>("[5,5.1]x[4,4.1]"));
  HybridSet<R> initial_set;
  initial_set.new_location(state1,initial_box);
  initial_set.new_location(state2, EmptySet<R>(2));
  initial_set.new_location(state3, EmptySet<R>(2));
  initial_set.new_location(state4, EmptySet<R>(2));
  cout << "initial_set.discrete_locations()=" << initial_set.locations() << endl;
  cout << "initial_set=" << initial_set << endl;
  
//  HybridSet<R> bounding_set;
//  bounding_set.new_location(state1,bounding_box);
//  bounding_set.new_location(state2,bounding_box);
//  bounding_set.new_location(state3,bounding_box);
//  bounding_set.new_location(state4,bounding_box);
//  cout << "bounding_set.discrete_locations()=" << bounding_set.locations() << endl;
//  cout << "bounding_set=" << bounding_set << endl;
//  cout << endl;

  cout << "Computing chainreachable set..." << endl;

  HybridSet<R> chainreach=hybrid_evolver.chainreach(automaton,initial_set);
  cout << "Reached " << chainreach
       << endl << endl;

  epsfstream eps; 

  eps.open("two_tank1-xv.eps",domain);
  //  eps << fill_colour(cyan) << activation;
  eps << line_style(true) << fill_colour(yellow) << chainreach[state1];
  eps << fill_colour(blue) << initial_set[state1];
  eps.close();

  eps.open("two_tank2-xv.eps",domain);
  //  eps << fill_colour(cyan) << activation;
  eps << line_style(true) << fill_colour(yellow) << chainreach[state2];
  eps << fill_colour(blue) << initial_set[state2];
  eps.close();

  
  eps.open("two_tank1.eps",domain,PlanarProjectionMap(3,2,0));
  //  eps << fill_colour(cyan) << activation;
  eps << line_style(true) << fill_colour(yellow) << chainreach[state1];
  eps << fill_colour(blue) << initial_set[state1];
  eps.close();
  
  return 0;
}
