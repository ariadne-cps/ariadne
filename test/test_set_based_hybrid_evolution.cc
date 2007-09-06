/***************************************************************************
 *            test_set_based_hybrid_evolution.cc
 *
 *  Copyright  2006  Alberto Casagrande,  Pieter Collins
 *  Email  casagrande@dimi.uniud.it  Pieter.Collins@cwi.nl
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

#include <iostream>
#include <fstream>
#include <string>

#include "ariadne.h"
#include "test_float.h"
#include "geometry/set_interface.h"
#include "geometry/set_reference.h"
#include "geometry/hybrid_set.h"
#include "geometry/rectangle.h"
#include "geometry/empty_set.h"
#include "geometry/polyhedral_set.h"
#include "system/affine_map.h"
#include "system/affine_vector_field.h"
#include "system/hybrid_automaton.h"
#include "system/set_based_hybrid_automaton.h"
#include "evaluation/applicator.h"
#include "evaluation/lohner_integrator.h"
#include "evaluation/affine_integrator.h"
#include "evaluation/hybrid_evolver.h"
#include "evaluation/set_based_hybrid_evolver.h"
#include "output/epsstream.h"
#include "output/logging.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Combinatoric;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;
using namespace Ariadne::Output;
using namespace std;

template<class R> int test_set_based_hybrid_evolution();
  
namespace Ariadne { namespace Evaluation { extern int verbosity; } }

template<class R>
int test_set_based_hybrid_evolution() 
{  
  //set_hybrid_evolver_verbosity(4);

  PolyhedralSet<R> space(Rectangle<R>("[-7.5,7.5]x[-7.5,7.5]"));
  AffineMap<R> identity(Matrix<R>("[1,0;0,-1]"),Vector<R>("[0,0]"));

  PolyhedralSet<R> invariant1(space);
  PolyhedralSet<R> invariant2(space);
  AffineVectorField<R> dynamic1(Matrix<R>("[-2,-1;1,-2]"),Vector<R>("[-1,0]"));
  AffineVectorField<R> dynamic2(Matrix<R>("[-2,-1; 1,-2]"),Vector<R>("[1,0]"));
  //AffineVectorField<R> dynamic2(Matrix<R>("[0,0;0,0]"),Vector<R>("[-1,-2]"));
  PolyhedralSet<R> activation11(Rectangle<R>("[-7.5,7.5]x[-3,-2]"));
  PolyhedralSet<R> activation21(Rectangle<R>("[-7.5,7.5]x[-7,-6]"));
  PolyhedralSet<R> activation12(Rectangle<R>("[-7.5,7.5]x[6,7]"));
  AffineMap<R> reset11(Matrix<R>("[1,0;0,-1]"),Vector<R>("[0,0]"));
  AffineMap<R> reset21(Matrix<R>("[1,0;0,-1]"),Vector<R>("[0,0]"));
  AffineMap<R> reset12(Matrix<R>("[1,0;0,-1]"),Vector<R>("[0,0]"));
   
  
  SetBasedHybridAutomaton<R> automaton("Affine automaton");
  id_type mode1_id=2;
  id_type mode2_id=3;
  const SetBasedDiscreteMode<R>& mode1=automaton.new_mode(mode1_id,dynamic1,invariant1);
  const SetBasedDiscreteMode<R>& mode2=automaton.new_mode(mode2_id,dynamic2,invariant2);
  id_type event1_id=5;
  id_type event2_id=7;
  const SetBasedDiscreteTransition<R>& transition11=automaton.new_transition(event1_id,mode1_id,mode1_id,reset11,activation11);
  const SetBasedDiscreteTransition<R>& transition21=automaton.new_transition(event1_id,mode2_id,mode1_id,reset21,activation21);
  const SetBasedDiscreteTransition<R>& transition12=automaton.new_transition(event2_id,mode1_id,mode2_id,reset12,activation12);
  
  cout << mode1 << " " << mode2 << endl;
  cout << transition11 << " " << transition21 << " " << transition12 << " " << endl;
  cout << endl;

  time_type maximum_step_size=0.125;
  time_type lock_to_grid_time=0.25;
  R maximum_set_radius=0.5;
  R grid_size=0.125;

  Applicator<R> apply(maximum_set_radius,grid_size);
  AffineIntegrator<R> integrator(maximum_step_size,lock_to_grid_time,maximum_set_radius); 
  SetBasedHybridEvolver<R> hybrid_evolver(apply,integrator);
  
  Grid<R> grid(Vector<R>("[0.25,0.25]"));
  FiniteGrid<R> finite_grid(grid,LatticeBlock("[-32,32]x[-32,32]"));
  Rectangle<R> bounding_box=finite_grid.extent();
  
  Rectangle<R> initial_rectangle1("[-6.96875,-6.9375]x[-6.96875,-6.9375]");
  Rectangle<R> initial_rectangle2("[6.9375,6.96875]x[6.9375,6.96875]");
  HybridGridMaskSet<R> initial_set;
  initial_set.new_location(mode1_id,finite_grid);
  initial_set.new_location(mode2_id,finite_grid);
  initial_set[mode1_id].adjoin_over_approximation(initial_rectangle1);
  //initial_set[mode2_id].adjoin_over_approximation(initial_rectangle2);
  cout << "initial_set=" << initial_set << endl;

  HybridGridMaskSet<R> bounding_set;
  bounding_set.new_location(mode1_id,finite_grid);
  bounding_set.new_location(mode2_id,finite_grid);
  bounding_set[mode1_id].adjoin_over_approximation(bounding_box);
  bounding_set[mode2_id].adjoin_over_approximation(bounding_box);
  cout << "bounding_set=" << bounding_set << endl;
  cout << endl;

  assert((bool)(subset(initial_set[mode1_id],bounding_set[mode1_id])));

  cout << "Computing continuous chainreach set" << endl;
  HybridGridMaskSet<R> continuous_chainreach=hybrid_evolver.continuous_chainreach(automaton,initial_set,bounding_set);
  cout << "Reached (" << continuous_chainreach[mode1_id].size() << "," << continuous_chainreach[mode2_id].size() << ") cells "
       << "out of (" << continuous_chainreach[mode1_id].capacity() << "," << continuous_chainreach[mode1_id].capacity() << ") "
       << " by continuous evolution" << endl << endl;
  
  epsfstream eps;
  eps.open("test_hybrid_evolution-1.eps",bounding_box.expand(0.5));
  eps << fill_colour(white) << bounding_box;
  eps << fill_colour(cyan) << automaton.mode(mode1_id).invariant();
  //eps << invariant1;
  eps << fill_colour(red) << continuous_chainreach[mode2_id];
  eps << fill_colour(green) << continuous_chainreach[mode1_id];
  eps << fill_colour(blue) << initial_set[mode1_id] << initial_set[mode2_id];
  eps.close();



  HybridGridMaskSet<R> initial_activated_set;
  initial_activated_set.new_location(mode1_id,finite_grid);
  initial_activated_set.new_location(mode2_id,finite_grid);
  cout << "empty initial_activated_set=" << initial_activated_set << endl;
  initial_activated_set[mode1_id].adjoin(over_approximation(Rectangle<R>("[-6,-5]x[-3.5,-1.0]"),grid));
  cout << "initial_activated_set=" << initial_activated_set << endl;
  cout << "initial_activated_set[mode1_id].size()=" << initial_activated_set[mode1_id].size() << " cells" << endl;
  cout << "Computing single discrete step" << endl;
  HybridGridMaskSet<R> discrete_reach=hybrid_evolver.discrete_step(automaton,initial_activated_set);
  cout << "Reached " << discrete_reach[mode2_id].size() << " cells by discrete step" << endl << endl;
  
  eps.open("test_hybrid_evolution-2.eps",bounding_box.expand(0.5));
  eps << fill_colour(white) << bounding_box;
  eps << fill_colour(cyan);
  eps << static_cast<const Polyhedron<R>&>(activation11);
  eps << activation11;
  eps << static_cast<const Polyhedron<R>&>(activation21);
  eps << activation21;
  eps << fill_colour(magenta);
  eps << static_cast<const Polyhedron<R>&>(activation12);
  eps << activation12;
  eps << fill_colour(red);
  eps << discrete_reach[mode2_id];
  eps << fill_colour(green);
  eps << discrete_reach[mode1_id];
  eps << fill_colour(blue);
  eps << initial_activated_set[mode1_id];
  eps.close();
  
  

  cout << "Computing chain reachable set" << endl;
  HybridGridMaskSet<R> chainreach=hybrid_evolver.chainreach(automaton,initial_set,bounding_set);
  cout << "Reached (" << chainreach[mode1_id].size() << "," << chainreach[mode2_id].size() << ") cells "
       << "out of (" << chainreach[mode1_id].capacity() << "," << chainreach[mode1_id].capacity() << ") "
       << endl << endl;
  
  eps.open("test_hybrid_evolution-3.eps",bounding_box.expand(0.5));
  eps << fill_colour(white) << bounding_box;
  eps << fill_colour(cyan) << activation11 << activation21;
  eps << fill_colour(magenta) << activation12;
  eps << line_style(true);
  eps << fill_colour(red) << chainreach[mode2_id];
  eps << fill_colour(yellow) << chainreach[mode1_id];
  eps << fill_colour(blue) << initial_set[mode1_id];
  eps.close();
  
  {
  cout << "\n\nUsing HybridSet arguments" << endl;
  HybridSet<R> initial_set;
  initial_set.new_location(mode1_id,initial_rectangle1);
  initial_set.new_location(mode2_id,EmptySet<R>(2));
  cout << "initial_set.discrete_locations()=" << initial_set.locations() << endl;
  cout << "initial_set[mode1_id]=" << initial_set[mode1_id] << endl;
  cout << "initial_set[mode21_id]=" << initial_set[mode2_id] << endl;

  HybridSet<R> bounding_set;
  bounding_set.new_location(mode1_id,bounding_box);
  bounding_set.new_location(mode2_id,bounding_box);
  cout << "bounding_set.discrete_locations()=" << bounding_set.locations() << endl;
  cout << "bounding_set[mode1_id]=" << bounding_set[mode1_id] << endl;
  cout << "bounding_set[mode21_id]=" << bounding_set[mode2_id] << endl;
  cout << endl;

  }
  
  return 0;
}


int main() {
  test_set_based_hybrid_evolution<Float>();
  cerr << "INCOMPLETE ";
  return 0;
}

