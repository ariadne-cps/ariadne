/***************************************************************************
 *            bouncing-ball-3var.cc
 *
 ****************************************************************************/


#include <iostream>
#include <fstream>
#include <string>

#include "ariadne.h"
#include "numeric/float64.h"
#include "geometry/set_interface.h"
#include "geometry/set_reference.h"
#include "geometry/hybrid_set.h"
#include "geometry/box.h"
#include "geometry/polyhedron.h"
#include "geometry/polyhedral_set.h"
#include "system/affine_map.h"
#include "system/affine_vector_field.h"
#include "system/hybrid_automaton.h"
#include "evaluation/evolution_parameters.h"
#include "evaluation/standard_applicator.h"
#include "evaluation/map_evolver.h"
#include "evaluation/vector_field_evolver.h"
#include "evaluation/lohner_integrator.h"
#include "evaluation/affine_integrator.h"
#include "evaluation/set_based_hybrid_evolver.h"
#include "output/epsstream.h"
#include "output/logging.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;
using namespace Ariadne::Output;
using namespace Ariadne::Combinatoric;
using namespace std;

template<class R> int bouncing_ball_automaton();
  
int main() {
  return bouncing_ball_automaton<Float64>();
}

template<class R>
int bouncing_ball_automaton() 
{
  typedef Zonotope<R> BS;
  
  set_hybrid_evolver_verbosity(4);

  Box<R> domain("[0,15]x[-20,20]x[0,10]");
  cout << "domain=" << domain << endl;

  AffineVectorField<R> dynamic(Matrix<R>("[0,1,0;0,0,0;0,0,0]"),Vector<R>("[0,-9.8,1.0]"));
  cout << "dynamic=" << dynamic << endl;
  RectangularSet<R> invariant(domain);
  cout << "invariant=" << invariant << endl;

  RectangularSet<R> activation(Box<R>("[0.0,0.01]x[-20,20]x[0,10]"));
  cout << "activation=" << activation << endl;

  AffineMap<R> reset(Matrix<R>("[1,0,0;0,-1,0;0,0,1]"),Vector<R>("[0,0,0]"));
  cout << "reset=" << reset << endl;
  
  HybridAutomaton<R> automaton("Bouncing ball automaton");
  DiscreteState mode1_id(0);
  const DiscreteMode<R>& mode1=automaton.new_mode(mode1_id,dynamic,invariant);
  DiscreteEvent event_id(5);
  const DiscreteTransition<R>& transition=automaton.new_transition(event_id,mode1_id,mode1_id,reset,activation);
  
  cout << mode1  <<  "\n" << transition << endl;

  time_type maximum_step_size=0.125;
  time_type lock_to_grid_time=0.5;
  R maximum_set_radius=0.25;
  R grid_size=0.125;

  EvolutionParameters<R> parameters;
  parameters.set_maximum_step_size(0.125);
  parameters.set_lock_to_grid_time(0.5);
  parameters.set_maximum_basic_set_radius(0.25);
  parameters.set_grid_length(0.125);

  StandardApplicator<R> applicator;
	LohnerIntegrator<R> lohner;
  SetBasedHybridEvolver<BS> hybrid_evolver(parameters,applicator,lohner);
  
  Grid<R> grid(Vector<R>("[0.25,0.25,1.0]"));
  Box<R> bounding_box(domain.neighbourhood(1));
  FiniteGrid<R> finite_grid(grid,bounding_box);
  
  Box<R> initial_box("[10,10.1]x[0,0.1]x[0,0.1]");
  HybridGridMaskSet<R> initial_set;
  initial_set.new_location(mode1_id,finite_grid);
  initial_set[mode1_id].adjoin_over_approximation(initial_box);
  cout << "initial_set.discrete_locations()=" << initial_set.locations() << endl;
  cout << "initial_set[mode1_id].size()=" << initial_set[mode1_id].size() << " cells out of " << initial_set[mode1_id].capacity() << endl;
  
  cout << "Computing chainreachable set..." << endl;

  HybridGridMaskSet<R> chainreach=hybrid_evolver.upper_reach(automaton,initial_set,2.0);
  cout << "Reached (" << chainreach[mode1_id].size() << ") cells "
       << "out of (" << chainreach[mode1_id].capacity() << ") "
       << endl << endl;

  epsfstream eps; 

  eps.open("bouncing_ball-3var-xv.eps",bounding_box);
  eps << fill_colour(white) << bounding_box;
  eps << fill_colour(cyan) << activation;
  eps << line_style(true) << fill_colour(yellow) << chainreach[mode1_id];
  eps << fill_colour(blue) << initial_set[mode1_id];
  eps.close();
  
  eps.open("bouncing_ball-3var-tx.eps",bounding_box,PlanarProjectionMap(3,2,0));
  eps << fill_colour(white) << bounding_box;
  eps << fill_colour(cyan) << activation;
  eps << line_style(true) << fill_colour(yellow) << chainreach[mode1_id];
  eps << fill_colour(blue) << initial_set[mode1_id];
  eps.close();
  
  return 0;
}
