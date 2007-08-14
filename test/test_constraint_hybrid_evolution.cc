/***************************************************************************
 *            test_constraint_hybrid_evolution.cc
 *
 *  Copyright  2007  Pieter Collins
 *  Email  Pieter.Collins@cwi.nl
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
#include "system/constraint_hybrid_automaton.h"
#include "evaluation/applicator.h"
#include "evaluation/lohner_integrator.h"
#include "evaluation/affine_integrator.h"
#include "evaluation/hybrid_evolver.h"
#include "evaluation/constraint_hybrid_evolver.h"
#include "evaluation/constraint_hybrid_evolver_plugin.h"
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

static const int verbosity = 0;

static const id_type mode1_id = 1;
static const id_type mode2_id = 2;
static const id_type mode3_id = 3;
static const id_type event20_id = 6;
static const id_type event12_id = 4;
static const id_type event23_id = 5;
  

template<class R>
ConstraintHybridEvolver<R> construct_evolver() 
{
  time_type maximum_step_size=0.125;
  time_type lock_to_grid_time=0.5;
  R maximum_set_radius=0.25;
  R grid_size=0.125;
  
  Applicator<R> apply(maximum_set_radius,grid_size);
  LohnerIntegrator<R> lohner_integrator(maximum_step_size,lock_to_grid_time,maximum_set_radius); 
  return ConstraintHybridEvolver<R>(apply,lohner_integrator);
}



template<class R>
ConstraintHybridAutomaton<R> construct_automaton() 
{
    AffineVectorField<R> dynamic1(Matrix<R>("[-1,-0.2;0.2,-1]"),Vector<R>("[5,0]"));
    AffineVectorField<R> dynamic2(Matrix<R>("[0,0; 0,0]"),Vector<R>("[0.75,1]"));
    AffineVectorField<R> dynamic3(Matrix<R>("[0,0; 0,0]"),Vector<R>("[1,-0.75]"));

    AffineMap<R> reset12(Matrix<R>("[0,1;-1,0]"),Vector<R>("[-2,3]"));
    LinearConstraint<R> guard12(Vector<R>("[1,0]"),Geometry::less,R(2));
    AffineMap<R> reset23(Matrix<R>("[1,0;0,1]"),Vector<R>("[0,-3]"));
    LinearConstraint<R> activation23(Vector<R>("[0,1]"),Geometry::less,R(2));
    LinearConstraint<R> invariant2(Vector<R>("[0,1]"),Geometry::less,R(2.5));
    
    ConstraintHybridAutomaton<R> automaton("");

    automaton.new_mode(mode1_id,dynamic1);
    automaton.new_mode(mode2_id,dynamic2);
    automaton.new_mode(mode3_id,dynamic3);
    automaton.new_invariant(event20_id,mode2_id,invariant2);
    automaton.new_forced_transition(event12_id,mode1_id,mode2_id,reset12,guard12);
    automaton.new_unforced_transition(event23_id,mode2_id,mode3_id,reset23,activation23);

    cout << "automaton = " << flush;
    cout << automaton << endl << endl;

    return automaton;
}





template<class R>
class TestConstraintHybridEvolution
{
 public:
  typedef Interval<R> I;

  ConstraintHybridAutomaton<R> automaton;
  ConstraintHybridEvolver<R> evolver;
  Polyhedron<R> guard;
  Polyhedron<R> activation;
  Polyhedron<R> invariant;
    
  TestConstraintHybridEvolution()
    : automaton(construct_automaton<R>()),
      evolver(construct_evolver<R>())
  {
  }


 int test_evolution() {
    time_type tr=2.0;
    time_type t=3.0;
    size_type n=2;

    id_type initial_discrete_mode = mode1_id;
    Zonotope<I,I> initial_basic_set(Point<R>("(0.5,0)"),Matrix<R>("[0.05,0.0; 0.0,0.05]"));
    HybridListSet< Zonotope<I,I> > initial_set(automaton.locations());
    initial_set.adjoin(initial_discrete_mode,initial_basic_set);
  

    
    cout << endl;

    Rectangle<R> bounding_box("[-3,3]x[-3,3]");
    Polyhedron<R> bounding_polyhedron(bounding_box);
    Polyhedron<R> guard_polyhedron("[[-1,0;-2]]");
    Polyhedron<R> activation_polyhedron("[[0,-1;-2]]");
    
    HybridListSet< Zonotope<I,I> > evolved_set=initial_set;
    evolved_set.clear();
    HybridListSet< Zonotope<I,I> > reached_set=initial_set;
    reached_set.clear();
    try {
      evolved_set=evolver.upper_evolve(automaton,initial_set,t,n);
      //reached_set=evolver.upper_reach(automaton,initial_set,tr,n);
    } 
    catch(const std::exception& e) {
      cout << "\nCaught: " << e.what() << endl;
    }
    //cout << "trace=" << evolver.trace() << endl << endl;

    cout << "initial_set = " << initial_set << endl;
    cout << "evolved_set = " << evolved_set << endl;
    cout << "reached_set = " << reached_set << endl;
    
    epsfstream eps;
    eps.open("test_constraint_hybrid_evolution-1.eps",bounding_box);
    eps << fill_colour(magenta) << closed_intersection(bounding_polyhedron,guard_polyhedron);
    eps << fill_colour(cyan) << closed_intersection(bounding_polyhedron,activation_polyhedron);

    for(uint i=0; i!=evolver.trace().size(); ++i) {
      switch(evolver.trace()[i].discrete_state()) {
        case mode1_id: eps << fill_colour(Output::green); break;
        case mode2_id: eps << fill_colour(Output::red); break;
        case mode3_id: eps << fill_colour(Output::blue); break;
      }
      eps << evolver.trace()[i].continuous_state_set();
    }      

    eps << fill_colour(Output::yellow) << initial_set[mode1_id];
    eps << reached_set[mode1_id];
    eps << evolved_set[mode1_id];
    eps << evolved_set[mode2_id];
    eps << evolved_set[mode3_id];

    
    /*
    Zonotope<I,I> iz(Point<I>("([-1,1],[-0.1,0.1])"),Matrix<I>("[1.0,0.2;1.0,-0.2]"));
    ListSet< Zonotope<I,I> > izls(iz);

    eps << iz;
    eps.set_fill_colour("red");
    eps << izls;
    */    


    eps.close();
    
    cout << automaton;

    return 0;
  }
  

  int test() {
    test_evolution();
    return 0;
  }

};



int main() {
  set_hybrid_evolver_verbosity(::verbosity);

  TestConstraintHybridEvolution<Float>().test();
  cerr << "INCOMPLETE ";
  return 0;
}

