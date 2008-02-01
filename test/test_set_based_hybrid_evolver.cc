/***************************************************************************
 *            test_set_based_hybrid_evolver.cc
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
#include "geometry/zonotope.h"
#include "geometry/empty_set.h"
#include "geometry/polyhedral_set.h"
#include "system/affine_map.h"
#include "system/affine_vector_field.h"
#include "system/hybrid_automaton.h"
#include "evaluation/evolution_parameters.h"
#include "evaluation/standard_applicator.h"
#include "evaluation/lohner_integrator.h"
#include "evaluation/affine_integrator.h"
#include "evaluation/map_evolver.h"
#include "evaluation/vector_field_evolver.h"
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

template<class R> int test_set_based_hybrid_evolver();
  
namespace Ariadne { namespace Evaluation { extern int verbosity; } }

template<class R>
class TestSetBasedHybridEvolver 
{  
  typedef Zonotope<R> BS;
  SetBasedHybridEvolver<BS>* hybrid_evolver;
  HybridAutomaton<R> automaton;
  HybridGridMaskSet<R> initial_set;
  HybridGridMaskSet<R> bounding_set;
  Box<R> bounding_box;
  DiscreteState mode1_id, mode2_id;

 public:
  TestSetBasedHybridEvolver() 
    : automaton(""),
      mode1_id(2),
      mode2_id(3)
  {
    hybrid_evolver_verbosity = 0;
    integrator_verbosity = 0;
    
    
    PolyhedralSet<R> space(Box<R>("[-7.5,7.5]x[-7.5,7.5]"));
    AffineMap<R> identity(Matrix<R>("[1,0;0,-1]"),Vector<R>("[0,0]"));
    
    PolyhedralSet<R> invariant1(space);
    PolyhedralSet<R> invariant2(space);
    AffineVectorField<R> dynamic1(Matrix<R>("[-2,-1;1,-2]"),Vector<R>("[-1,0]"));
    AffineVectorField<R> dynamic2(Matrix<R>("[-2,-1; 1,-2]"),Vector<R>("[1,0]"));
    //AffineVectorField<R> dynamic2(Matrix<R>("[0,0;0,0]"),Vector<R>("[-1,-2]"));
    PolyhedralSet<R> activation11(Box<R>("[-7.5,7.5]x[-3,-2]"));
    PolyhedralSet<R> activation21(Box<R>("[-7.5,7.5]x[-7,-6]"));
    PolyhedralSet<R> activation12(Box<R>("[-7.5,7.5]x[6,7]"));
    AffineMap<R> reset11(Matrix<R>("[1,0;0,-1]"),Vector<R>("[0,0]"));
    AffineMap<R> reset21(Matrix<R>("[1,0;0,-1]"),Vector<R>("[0,0]"));
    AffineMap<R> reset12(Matrix<R>("[1,0;0,-1]"),Vector<R>("[0,0]"));
    
    
    automaton=HybridAutomaton<R>("Affine automaton");
    mode1_id=DiscreteState(2);
    mode2_id=DiscreteState(3);
    const DiscreteMode<R>& mode1=automaton.new_mode(mode1_id,dynamic1,invariant1);
    const DiscreteMode<R>& mode2=automaton.new_mode(mode2_id,dynamic2,invariant2);
    DiscreteEvent event1_id(5);
    DiscreteEvent event2_id(7);
    const DiscreteTransition<R>& transition11=automaton.new_transition(event1_id,mode1_id,mode1_id,reset11,activation11);
    const DiscreteTransition<R>& transition21=automaton.new_transition(event1_id,mode2_id,mode1_id,reset21,activation21);
    const DiscreteTransition<R>& transition12=automaton.new_transition(event2_id,mode1_id,mode2_id,reset12,activation12);
    
    cout << mode1 << " " << mode2 << endl;
    cout << transition11 << " " << transition21 << " " << transition12 << " " << endl;
    cout << endl;
    
    time_type maximum_step_size=0.125;
    time_type lock_to_grid_time=0.25;
    R maximum_basic_set_radius=0.5;
    R grid_length=0.25;
    
    EvolutionParameters<R> parameters;
    parameters.set_maximum_basic_set_radius(maximum_basic_set_radius);
    parameters.set_grid_length(grid_length);
    parameters.set_maximum_step_size(maximum_step_size);
    parameters.set_lock_to_grid_time(lock_to_grid_time);
    parameters.set_verbosity(0);
    
    StandardApplicator<R> applicator;
    ApplicatorInterface<BS>& applicator_interface=applicator;
    AffineIntegrator<R> affine_integrator;
    IntegratorInterface<BS>& integrator_interface=affine_integrator;
    MapEvolver<BS> discrete_time_evolver(parameters,applicator_interface);
    VectorFieldEvolver<BS> continuous_time_evolver(parameters,integrator_interface);
    hybrid_evolver=new SetBasedHybridEvolver<BS> (parameters,applicator_interface,integrator_interface);
    
    Grid<R> grid(Vector<R>("[0.25,0.25]"));
    FiniteGrid<R> finite_grid(grid,LatticeBlock("[-32,32]x[-32,32]"));
    bounding_box=finite_grid.extent();
    RectangularSet<R> bounding_rectangle(bounding_box);
    
    RectangularSet<R> initial_rectangle1(Box<R>("[-6.96875,-6.9375]x[-6.96875,-6.9375]"));
    RectangularSet<R> initial_rectangle2(Box<R>("[6.9375,6.96875]x[6.9375,6.96875]"));

    initial_set.new_location(mode1_id,finite_grid);
    initial_set.new_location(mode2_id,finite_grid);
    initial_set[mode1_id].adjoin_outer_approximation(initial_rectangle1);
    //initial_set[mode2_id].adjoin_outer_approximation(initial_rectangle2);
    cout << "initial_set=" << initial_set << endl;
    
    bounding_set.new_location(mode1_id,finite_grid);
    bounding_set.new_location(mode2_id,finite_grid);
    bounding_set[mode1_id].adjoin_over_approximation(bounding_box);
    bounding_set[mode2_id].adjoin_over_approximation(bounding_box);
    cout << "bounding_set=" << bounding_set << endl;
    cout << endl;
    
    assert((bool)(subset(initial_set[mode1_id],bounding_set[mode1_id])));
    
  }

  void test_upper_reach() {  
    cout << "Computing timed reachable set" << endl;
    HybridGridMaskSet<R> upper_reach=hybrid_evolver->upper_reach(automaton,initial_set,Rational(1));
    cout << "Reached (" << upper_reach[mode1_id].size() << "," << upper_reach[mode2_id].size() << ") cells "
         << "out of (" << upper_reach[mode1_id].capacity() << "," << upper_reach[mode1_id].capacity() << ") "
         << endl << endl;
  }

  void test_chain_reach() {  
    cout << "Computing chain reachable set" << endl;
    HybridGridMaskSet<R> chainreach=hybrid_evolver->chainreach(automaton,initial_set);
    cout << "Reached (" << chainreach[mode1_id].size() << "," << chainreach[mode2_id].size() << ") cells "
         << "out of (" << chainreach[mode1_id].capacity() << "," << chainreach[mode1_id].capacity() << ") "
         << endl << endl;

    epsfstream eps;
    eps.open("test_hybrid_evolver-chain_reach.eps",bounding_box.neighbourhood(0.5));
    eps << fill_colour(white) << bounding_box;
    //eps << fill_colour(cyan) << activation11 << activation21;
    //eps << fill_colour(magenta) << activation12;
    eps << line_style(true);
    eps << fill_colour(red) << chainreach[mode2_id];
    eps << fill_colour(yellow) << chainreach[mode1_id];
    eps << fill_colour(blue) << initial_set[mode1_id];
    eps.close();
  }
  
  void test() {
    ARIADNE_TEST_CALL(test_upper_reach());
    //ARIADNE_TEST_CALL(test_chain_reach());
  }

};


int main(int nargs, const char* args[]) 
{
  int hybrid_evolver_verbosity = 0;
  int integrator_verbosity = 0;
  if(nargs>1) {
    hybrid_evolver_verbosity = std::atoi(args[1]);
  }
  if(nargs>2) {
    integrator_verbosity = std::atoi(args[2]);
  }
  set_hybrid_evolver_verbosity(hybrid_evolver_verbosity);
  set_integrator_verbosity(integrator_verbosity);
  TestSetBasedHybridEvolver<Flt>().test();
  cerr << "INCOMPLETE ";
  return ARIADNE_TEST_FAILURES;
}

