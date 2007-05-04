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
#include "geometry/rectangle.h"
#include "geometry/polyhedron.h"
#include "geometry/polyhedral_set.h"
#include "system/affine_map.h"
#include "system/affine_vector_field.h"
#include "system/hybrid_automaton.h"
#include "evaluation/applicator.h"
#include "evaluation/lohner_integrator.h"
#include "evaluation/affine_integrator.h"
#include "evaluation/hybrid_evolver.h"
#include "output/epsfstream.h"
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
  
  Rectangle<R> r("[0,15]x[-20,20]x[0,100]");
  cout << "r=" << r << endl;

  AffineVectorField<R> dynamic(Matrix<R>("[0,-1,0;0,0,0;0,0,0]"),Vector<R>("[0,-9.8,1.0]"));
  cout << "dynamic=" << dynamic << endl;
  PolyhedralSet<R> invariant(r);
  cout << "invariant=" << invariant << endl;

  PolyhedralSet<R> activation(Rectangle<R>("[0.0,0.01]x[-20,20]x[0,100]"));
  cout << "activation=" << activation << endl;

  AffineMap<R> reset(Matrix<R>("[1,0,0;0,-1,0;0,0,1]"),Vector<R>("[0,0,0]"));
  cout << "reset=" << reset << endl;
  
  HybridAutomaton<R> automaton("Bouncing ball automaton");
  id_type mode1_id=0;
  const DiscreteMode<R>& mode1=automaton.new_mode(mode1_id,dynamic,invariant);
  id_type event_id=5;
  const DiscreteTransition<R>& transition=automaton.new_transition(event_id,mode1_id,mode1_id,reset,activation);
  
  cout << mode1  <<  "\n" << transition << endl;
  cout << automaton.invariant() << endl;

	time_type maximum_step_size=0.125;
  time_type lock_to_grid_time=0.5;
  R maximum_set_radius=0.25;

  Applicator<R> apply;
  AffineIntegrator<R> integrator(maximum_step_size,lock_to_grid_time,maximum_set_radius); 
  HybridEvolver<R> hybrid_evolver(apply,integrator);
  
  Grid<R> grid(Vector<R>("[0.25,0.25,1.0]"));
  FiniteGrid<R> finite_grid(grid,LatticeBlock("[-10,100]x[-100,100]x[-5,105]"));
  
  Rectangle<R> initial_rectangle("[10,10.1]x[0,0.1]x[0,0.1]");
  HybridGridMaskSet<R> initial_set;
  initial_set.new_location(mode1_id,finite_grid);
  initial_set[mode1_id].adjoin_over_approximation(initial_rectangle);
  cout << "initial_set.discrete_locations()=" << initial_set.locations() << endl;
  cout << "initial_set[mode1_id].size()=" << initial_set[mode1_id].size() << " cells out of " << initial_set[mode1_id].capacity() << endl;
  
  Rectangle<R> bounding_box("[-1,16]x[-21,21]x[-1,101]");
  HybridGridMaskSet<R> bounding_set;
  bounding_set.new_location(mode1_id,finite_grid);
  bounding_set[mode1_id].adjoin_over_approximation(bounding_box);
  cout << "bounding_set.discrete_locations()=" << bounding_set.locations() << endl;
  cout << "bounding_set[mode1_id].size()=" << bounding_set[mode1_id].size() << " cells out of " << bounding_set[mode1_id].capacity() << endl;
  cout << endl;

  assert((bool)(subset(initial_set[mode1_id],bounding_set[mode1_id])));

	cout << "Computing chainreachable set..." << endl;

  HybridGridMaskSet<R> chainreach=hybrid_evolver.chainreach(automaton,initial_set,bounding_set);
  cout << "Reached (" << chainreach[mode1_id].size() << ") cells "
       << "out of (" << chainreach[mode1_id].capacity() << ") "
       << endl << endl;

	epsfstream eps; 
  eps.open("bouncing-ball-3var-cc.eps",bounding_box.expand(0.5));
  eps.set_fill_colour("white");
  eps << bounding_box;
  eps.set_fill_colour("cyan");
  eps << static_cast<const Polyhedron<R>&>(activation);
  eps.set_fill_colour("magenta");
  eps.set_line_style(true);
  eps.set_fill_colour("yellow");
  eps << chainreach[mode1_id];
  eps.set_fill_colour("blue");
  eps << initial_set[mode1_id];
  eps.close();
  
  return 0;
}
