/****************************************************************************
 * simple.cc
 * Simple example of hybrid automata with locations of different dimensions
 * (c) Davide Bresolin davide.bresolin@univr.it , 2008
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
#include "evaluation/evolution_parameters.h"
#include "evaluation/standard_applicator.h"
#include "evaluation/kuhn_applicator.h"
#include "evaluation/lohner_integrator.h"
#include "evaluation/affine_integrator.h"
#include "evaluation/set_based_hybrid_evolver.h"
#include "output/epsstream.h"
#include "output/textstream.h"
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
  typedef Zonotope<R> BS;

  Box<R> domain("[0,1]x[0,1]"); //Setup domain
  cout << "domain=" << domain << endl;

  double a = -0.02;
  double b = 0.3;
  double T = 4.0;
  double h = 5.5;
  double Delta = 0.05;
  double H = 8.0;
  
  double A[4]={a,b,0,0};
  double v[2]={0,1.0/T};
  AffineVectorField<R> dynamic1(Matrix<R>(2,2,A),Vector<R>(2,v));
  cout << "dynamic1=" << dynamic1 << endl;

  RectangularSet<R> inv1(Box<R>("[0,5.55]x[0,1]")); 
  cout << "invariant1=" << inv1 << endl;

  double B[1]={a};
  double u[1]={b};
  AffineVectorField<R> dynamic2(Matrix<R>(1,1,B),Vector<R>(1,u));
  cout << "dynamic2=" << dynamic2 << endl;

  RectangularSet<R> inv2(Box<R>("[0,5.55]")); 
  cout << "invariant2=" << inv2 << endl;
  
  HybridAutomaton<R> automaton("Simple Automata");
  // Building continues modes
  DiscreteState state1(1);
  DiscreteState state2(2);
  cout << "states={" << state1 << " ; " << state2 << "}" << endl;

  const DiscreteMode<R>& mode1=automaton.new_mode(state1,dynamic1,inv1);
  const DiscreteMode<R>& mode2=automaton.new_mode(state2,dynamic2,inv2);  
  cout << "modes={" << mode1 << " ; " << mode2 << "}" << endl;

  RectangularSet<R> guard12(Box<R>("[0,8]x[1,1.5]"));
  cout << "guard12=" << guard12 << endl;

  AffineMap<R> reset12(Matrix<R>("[1,0]"),Vector<R>("[0]"));
  cout << "reset12=" << reset12 << endl;

  RectangularSet<R> guard21(Box<R>("[0,0]"));
  cout << "guard21=" << guard21 << endl;

  AffineMap<R> reset21(Matrix<R>("[1;0]"),Vector<R>("[0,0]"));
  cout << "reset21=" << reset21 << endl;

  // Form the forced transitions 
  DiscreteEvent event12(5);
  DiscreteEvent event21(6);
  cout << "events={" << event12 << " " << event21 << "}" << endl;

  const DiscreteTransition<R>& transition12=automaton.new_transition(event12,state1,state2,reset12,guard12);
  cout << "transition1=" << transition12 << endl;

  const DiscreteTransition<R>& transition21=automaton.new_transition(event21,state2,state1,reset21,guard21);
  // cout << "transition21=" << transition21 << endl;

  EvolutionParameters<R> parameters;
  //cout << "parameters.minimum_step_size()=" << flush; cout << Rational(parameters.minimum_step_size()) << endl;
  parameters.set_maximum_step_size(0.125);
  parameters.set_lock_to_grid_time(0.5);
  parameters.set_maximum_basic_set_radius(0.25);
  parameters.set_grid_length(0.125);
  parameters.set_verbosity(9);
  set_geometry_verbosity(0);
  cout << "parameters=" << parameters << endl;

  // BUG: StandardApplicator does not support locations with different dimensions
  StandardApplicator<R> apply;
  // KuhnApplicator works:
  // KuhnApplicator<R> apply(3);
  AffineIntegrator<R> integrator;
  SetBasedHybridEvolver<BS> hybrid_evolver(parameters,apply,integrator);
  
  RectangularSet<R> bounding_box(domain.neighbourhood(1));
  cout << "bounding_box=" << bounding_box << endl;
  
  RectangularSet<R> initial_box(Box<R>("[0,0]x[0,0]"));
  HybridSet<R> initial_set;
  initial_set.new_location(state1,initial_box);
  initial_set.new_location(state2,EmptySet<R>(1));
  cout << "initial_set.discrete_locations()=" << initial_set.locations() << endl;
  cout << "initial_set=" << initial_set << endl;
  
  R time=16.0;
  cout << "Computing reachable set for "<< time <<"seconds..." << endl;

  HybridSet<R> reach=hybrid_evolver.upper_reach(automaton,initial_set,time);
  cout << "Reached set=" << reach
       << endl << endl;

  epsfstream eps; 

  eps.open("simple1.eps",domain);
  //  eps << fill_colour(cyan) << activation;
  eps << line_style(true) << fill_colour(yellow) << reach[state1];
  eps << fill_colour(blue) << initial_set[state1];
  eps.close();
    
  textfstream txt;
  
  txt.open("simple1.txt");
  txt << reach[state1];
  txt.close();
  
  txt.open("simple2.txt");
  txt << reach[state2];
  txt.close();
  
  return 0;
}
