/***************************************************************************
 *            test_hybrid_evolution.cc
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

template<class R> int test_hybrid_evolution();
  
namespace Ariadne { namespace Evaluation { extern int verbosity; } }

int main() {
  test_hybrid_evolution<Float>();
  cerr << "INCOMPLETE ";
  return 0;
}

template<class R>
int test_hybrid_evolution() 
{  
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
   
  
  HybridAutomaton<R> automaton("Affine automaton");
  id_type mode1_id=2;
  id_type mode2_id=3;
  const DiscreteMode<R>& mode1=automaton.new_mode(mode1_id,dynamic1,invariant1);
  const DiscreteMode<R>& mode2=automaton.new_mode(mode2_id,dynamic2,invariant2);
  id_type event1_id=5;
  id_type event2_id=7;
  const DiscreteTransition<R>& transition11=automaton.new_transition(event1_id,mode1_id,mode1_id,reset11,activation11);
  const DiscreteTransition<R>& transition21=automaton.new_transition(event1_id,mode2_id,mode1_id,reset21,activation21);
  const DiscreteTransition<R>& transition12=automaton.new_transition(event2_id,mode1_id,mode2_id,reset12,activation12);
  
  cout << mode1 << " " << mode2 << endl;
  cout << transition11 << " " << transition21 << " " << transition12 << " " << endl;
  cout << endl;

  time_type maximum_step_size=0.125;
  time_type lock_to_grid_time=0.25;
  R maximum_set_radius=0.5;

  Applicator<R> apply;
  AffineIntegrator<R> integrator(maximum_step_size,lock_to_grid_time,maximum_set_radius); 
  HybridEvolver<R> hybrid_evolver(apply,integrator);
  
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
  cout << "initial_set.discrete_locations()=" << initial_set.locations() << endl;
  cout << "initial_set[mode1_id].size()=" << initial_set[mode1_id].size() << " cells out of " << initial_set[mode1_id].capacity() << endl;
  cout << "initial_set[mode21_id].size()=" << initial_set[mode2_id].size() << " cells out of " << initial_set[mode2_id].capacity() << endl;

  HybridGridMaskSet<R> bounding_set;
  bounding_set.new_location(mode1_id,finite_grid);
  bounding_set.new_location(mode2_id,finite_grid);
  bounding_set[mode1_id].adjoin_over_approximation(bounding_box);
  bounding_set[mode2_id].adjoin_over_approximation(bounding_box);
  cout << "bounding_set.discrete_locations()=" << bounding_set.locations() << endl;
  cout << "bounding_set[mode1_id].size()=" << bounding_set[mode1_id].size() << " cells out of " << bounding_set[mode1_id].capacity() << endl;
  cout << "bounding_set[mode21_id].size()=" << bounding_set[mode2_id].size() << " cells out of " << bounding_set[mode2_id].capacity() << endl;
  cout << endl;

  assert((bool)(subset(initial_set[mode1_id],bounding_set[mode1_id])));

  cout << "Computing continuous chainreach set" << endl;
  HybridGridMaskSet<R> continuous_chainreach=hybrid_evolver.continuous_chainreach(automaton,initial_set,bounding_set);
  cout << "Reached (" << continuous_chainreach[mode1_id].size() << "," << continuous_chainreach[mode2_id].size() << ") cells "
       << "out of (" << continuous_chainreach[mode1_id].capacity() << "," << continuous_chainreach[mode1_id].capacity() << ") "
       << " by continuous evolution" << endl << endl;
  
  epsfstream eps;
  eps.open("test_hybrid_evolution-1.eps",bounding_box.expand(0.5));
  eps.set_fill_colour("white");
  eps << bounding_box;
  eps.set_fill_colour("red");
  eps << continuous_chainreach[mode2_id];
  eps.set_fill_colour("green");
  eps << continuous_chainreach[mode1_id];
  eps.set_fill_colour("blue");
  eps << initial_set[mode1_id];
  eps << initial_set[mode2_id];
  eps.close();


  HybridGridCellListSet<R> initial_activated_set;
  initial_activated_set.new_location(mode1_id,grid);
  initial_activated_set.new_location(mode2_id,grid);
  initial_activated_set[mode1_id].adjoin(over_approximation(Rectangle<R>("[-6,-5]x[-3.5,-1.0]"),grid));
  cout << "initial_activated_set[mode1_id].size()=" << initial_activated_set[mode1_id].size() << " cells" << endl;
  cout << "Computing single discrete step" << endl;
  HybridGridCellListSet<R> discrete_reach=hybrid_evolver.discrete_step(automaton,initial_activated_set);
  cout << "Reached " << discrete_reach[mode2_id].size() << " cells by discrete step" << endl << endl;
  
  eps.open("test_hybrid_evolution-2.eps",bounding_box.expand(0.5));
  eps.set_fill_colour("white");
  eps << bounding_box;
  eps.set_fill_colour("cyan");
  eps << static_cast<Polyhedron<R>&>(activation11);
  eps << static_cast<Polyhedron<R>&>(activation21);
  eps.set_fill_colour("magenta");
  eps << static_cast<Polyhedron<R>&>(activation12);
  eps.set_fill_colour("red");
  eps << discrete_reach[mode2_id];
  eps.set_fill_colour("green");
  eps << discrete_reach[mode1_id];
  eps.set_fill_colour("blue");
  eps << initial_activated_set[mode1_id];
  eps.close();
  
  

  HybridGridMaskSet<R> chainreach=hybrid_evolver.chainreach(automaton,initial_set,bounding_set);
  cout << "Reached (" << chainreach[mode1_id].size() << "," << chainreach[mode2_id].size() << ") cells "
       << "out of (" << chainreach[mode1_id].capacity() << "," << chainreach[mode1_id].capacity() << ") "
       << endl << endl;
  
  eps.open("test_hybrid_evolution-3.eps",bounding_box.expand(0.5));
  eps.set_fill_colour("white");
  eps << bounding_box;
  eps.set_fill_colour("cyan");
  eps << static_cast<Polyhedron<R>&>(activation11);
  eps << static_cast<Polyhedron<R>&>(activation21);
  eps.set_fill_colour("magenta");
  eps << static_cast<Polyhedron<R>&>(activation12);
  eps.set_line_style(true);
  eps.set_fill_colour("red");
  eps << chainreach[mode2_id];
  eps.set_fill_colour("yellow");
  eps << chainreach[mode1_id];
  eps.set_fill_colour("blue");
  eps << initial_set[mode1_id];
  eps.close();
  
  
  return 0;
}
