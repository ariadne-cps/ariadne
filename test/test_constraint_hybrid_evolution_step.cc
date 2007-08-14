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

static const int verbosity = 0;

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Combinatoric;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;
using namespace Ariadne::Output;
using namespace std;

static const id_type mode1_id = 1;
static const id_type mode2_id = 2;
static const id_type event3_id = 3;
static const id_type event4_id = 4;
static const id_type event5_id = 5;
static const id_type event6_id = 6;
  
template<class R>
ConstraintHybridEvolverPlugin<R> construct_evolver_plugin() 
{
  time_type maximum_step_size=0.125;
  time_type lock_to_grid_time=0.5;
  R maximum_set_radius=0.25;
  R grid_size=0.125;
  
  Applicator<R> apply(maximum_set_radius,grid_size);
  LohnerIntegrator<R> lohner_integrator(maximum_step_size,lock_to_grid_time,maximum_set_radius); 
  return ConstraintHybridEvolverPlugin<R>(apply,lohner_integrator);
}


template<class R>
ConstraintHybridAutomaton<R> construct_automaton() 
{
    AffineVectorField<R> dynamic1(Matrix<R>("[-0.3,-0.4;0.2,0.2]"),Vector<R>("[2,0.5]"));
    AffineVectorField<R> dynamic2(Matrix<R>("[0,0; 0,0]"),Vector<R>("[0.75,1]"));

    AffineMap<R> reset312(Matrix<R>("[0,1;-1,0]"),Vector<R>("[-2,3]"));
    LinearConstraint<R> guard312(Vector<R>("[1,0]"),Geometry::less,R(2));
    AffineMap<R> reset412(Matrix<R>("[0,1;-1,0]"),Vector<R>("[-2,1]"));
    LinearConstraint<R> guard412(Vector<R>("[0,-1]"),Geometry::less,R(2));
    AffineMap<R> reset512(Matrix<R>("[1,0;0,1]"),Vector<R>("[0,-3]"));
    LinearConstraint<R> activation512(Vector<R>("[0,1]"),Geometry::less,R(2));
    LinearConstraint<R> invariant61(Vector<R>("[0,1]"),Geometry::less,R(2.5));
    
    ConstraintHybridAutomaton<R> automaton("");

    automaton.new_mode(mode1_id,dynamic1);
    automaton.new_mode(mode2_id,dynamic2);
    automaton.new_forced_transition(event3_id,mode1_id,mode2_id,reset312,guard312);
    automaton.new_forced_transition(event4_id,mode1_id,mode2_id,reset412,guard412);
    automaton.new_unforced_transition(event5_id,mode1_id,mode2_id,reset512,activation512);
    automaton.new_invariant(event6_id,mode1_id,invariant61);

    cout << "automaton = " << flush;
    cout << automaton << endl << endl;

    return automaton;
}




template<class R>
class TestConstraintHybridEvolverPlugin
{
  typedef Interval<R> I;
  typedef TimeModelHybridBasicSet< Zonotope<I,I> > timed_set_type;
 public:
  ConstraintHybridAutomaton<R> automaton;
  ConstraintHybridEvolverPlugin<R> evolver;
  Rectangle<R> bounding_box;
    

  TestConstraintHybridEvolverPlugin()
    : automaton(construct_automaton<R>()), 
      evolver(construct_evolver_plugin<R>()),
      bounding_box("[-3,3]x[-3,3]")
  { }

  epsfstream& write_invariants(epsfstream& eps) const {
    Rectangle<R> guard3_rectangle("[2,3]x[-3,3]");
    Rectangle<R> guard4_rectangle("[-3,3]x[-3,-2]");
    Rectangle<R> activation5_rectangle("[-3,3]x[2,3]");
    Rectangle<R> invariant6_rectangle("[-3,3]x[2.5,3]");
    eps << fill_colour(cyan) << activation5_rectangle;
    eps << fill_colour(black) << invariant6_rectangle;
    eps << fill_colour(magenta) << guard4_rectangle;
    eps << fill_colour(magenta) << guard3_rectangle;
    return eps;
  }

  void plot(const char* name, timed_set_type initial_set, const std::vector<timed_set_type> reach_set) const
  {
    epsfstream eps;
    std::string filename=std::string("test_constraint_hybrid_evolution_step-")+name+".eps";
    eps.open(filename.c_str(),bounding_box);
    write_invariants(eps);
    assert(reach_set.size()>=1);
    eps << fill_colour(green);
    eps << reach_set.back().continuous_state_set();
    eps << fill_colour(white);
    for(size_type i=0; i!=reach_set.size()-1; ++i) {
      eps << reach_set[i].continuous_state_set();
    }
    eps << fill_colour(yellow);
    eps << initial_set.continuous_state_set();
    eps.close();
    return;
  }

  int test_reachability_step() {
    cout << "\n" << __FUNCTION__ << "\n";
    time_type step_size=0.25;
    timed_set_type initial_set(1,Zonotope<I,I>("{(0.5,-0.4),[0.1,0,0.01;0,0.1,0.01]}"));
    std::vector<timed_set_type> reach_set=evolver.evolution_step(automaton,initial_set,step_size,upper_semantics,compute_reachable_set);
    cout << "\nreach_set=" << reach_set << endl;
    plot("reach",initial_set,reach_set);
    return 0;
  }
  
  int test_guard_step() {
    cout << "\n" << __FUNCTION__ << "\n";
    time_type step_size=0.25;
    timed_set_type initial_set(1,Zonotope<I,I>("{(1.90,-0.4),[0.05,0,0.01;0,0.05,0.01]}"));
    std::vector<timed_set_type> reach_set=evolver.evolution_step(automaton,initial_set,step_size,upper_semantics,compute_reachable_set);
    cout << "\nreach_set=" << reach_set << endl;
    plot("guard",initial_set,reach_set);
    return 0;
  }
  
  int test_activation_step() {
    cout << "\n" << __FUNCTION__ << "\n";
    time_type step_size=0.25;
    timed_set_type initial_set(1,Zonotope<I,I>(Rectangle<R>("[-1.45,-1.35]x[2.05,2.15]")));
    std::vector<timed_set_type> reach_set=evolver.evolution_step(automaton,initial_set,step_size,upper_semantics,compute_reachable_set);
    cout << "\nreach_set=" << reach_set << endl;
    plot("activation",initial_set,reach_set);
    return 0;
  }
  
  int test_invariant_step() {
    cout << "\n" << __FUNCTION__ << "\n";
    time_type step_size=0.25;
    TimeModelHybridBasicSet< Zonotope<I,I> > initial_set(1,Zonotope<I,I>(Rectangle<R>("[1.35,1.45]x[2.35,2.45]")));
    std::vector<timed_set_type> reach_set=evolver.evolution_step(automaton,initial_set,step_size,upper_semantics,compute_reachable_set);
    cout << "\nreach_set=" << reach_set << endl;
    plot("invariant",initial_set,reach_set);
    return 0;
  }
  
  int test_initial_guard_step() {
    cout << "\n" << __FUNCTION__ << "\n";
    time_type step_size=0.25;
    TimeModelHybridBasicSet< Zonotope<I,I> > initial_set(1,Zonotope<I,I>(Rectangle<R>("[1.95,2.05]x[1.35,1.45]")));
    std::vector<timed_set_type> reach_set=evolver.evolution_step(automaton,initial_set,step_size,upper_semantics,compute_reachable_set);
    cout << "\nreach_set=" << reach_set << endl;
    plot("initial_guard",initial_set,reach_set);
    return 0;
  }
  
  int test_corner_collision_step() {
    cout << "\n" << __FUNCTION__ << "\n";
    time_type step_size=0.25;
    TimeModelHybridBasicSet< Zonotope<I,I> > initial_set(1,Zonotope<I,I>(Rectangle<R>("[1.65,1.75]x[2.35,2.45]")));
    std::vector<timed_set_type> subdivided_set=evolver.evolution_step(automaton,initial_set,step_size,upper_semantics,compute_reachable_set);
    std::vector<timed_set_type> reach_set;
    cout << "\nsubdivided_set=" << subdivided_set << endl;
    for(uint i=0; i!=subdivided_set.size(); ++i) { 
      std::vector<timed_set_type> reach_subset=evolver.evolution_step(automaton,subdivided_set[i],step_size,upper_semantics,compute_reachable_set);
      reach_set.insert(reach_set.end(),reach_subset.begin(),reach_subset.end());
    }
    cout << "reach_set=" << reach_set << endl;
    plot("corner_collision",initial_set,reach_set);
    return 0;
  }
  
  int test() {
    test_reachability_step();
    test_guard_step();
    test_activation_step();
    test_invariant_step();
    test_initial_guard_step();
    //test_corner_collision_step();
    return 0;
  }

};



int main() {
  set_hybrid_evolver_verbosity(::verbosity);

  TestConstraintHybridEvolverPlugin<Float>().test();
  cerr << "INCOMPLETE ";
  return 0;
}

