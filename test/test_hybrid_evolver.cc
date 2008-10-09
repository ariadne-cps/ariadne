/***************************************************************************
 *            test_hybrid_evolver.cc
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

#include <iostream>
#include <fstream>
#include <string>

#include "ariadne.h"
#include "test_float.h"
#include "geometry/hybrid_set.h"
#include "geometry/zonotope.h"
#include "geometry/empty_set.h"
#include "geometry/taylor_set.h"
#include "system/affine_map.h"
#include "system/affine_vector_field.h"
#include "system/hybrid_automaton.h"
#include "evaluation/evolution_parameters.h"
#include "evaluation/standard_applicator.h"
#include "evaluation/standard_bounder.h"
#include "evaluation/affine_integrator.h"
#include "evaluation/standard_satisfier.h"
#include "evaluation/standard_subdivider.h"
#include "evaluation/cascade_reducer.h"
#include "evaluation/hybrid_evolver.h"
#include "output/epsstream.h"
#include "output/logging.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

template<class R> int test_set_based_hybrid_evolver();
template<class R> int test_constraint_based_hybrid_evolver();
  
/*

template<class R>
class TestSetBasedHybridEvolver 
{  
  typedef HybridAutomaton<R> Sys;
  typedef Zonotope<R> ES;
  typedef HybridBasicSet< Zonotope<R> > HES;
  typedef ListSet< HybridBasicSet< Zonotope<R> > > HLS;

  Evolver<HybridAutomaton<R>, Zonotope<R> > hybrid_evolver;
  HybridAutomaton<R> automaton;
  HybridBasicSet< Zonotope<R> > initial_set;
  Box<R> bounding_box;
  DiscreteState mode1_id, mode2_id;

 public:
  TestHybridEvolver() 
    : automaton(""),
      initial_set(DiscreteState(2),Zonotope<R>(2)),
      bounding_box(Box<R>("[-4,4]x[-4,4]")),
      mode1_id(2),
      mode2_id(3)
  {
    hybrid_evolver_verbosity = 2;
    integrator_verbosity = 0;
    
    
    PolyhedralSet<R> spaceset(Box<R>("[-7.5,7.5]x[-7.5,7.5]"));
    AffineFunction<R> identity(Matrix<R>("[1,0;0,-1]"),Vector<R>("[0,0]"));
    
    PolyhedralSet<R> invariantset1(spaceset);
    PolyhedralSet<R> invariantset2(spaceset);
    AffineFunction<R> dynamic1(Matrix<R>("[-2,-1;1,-2]"),Vector<R>("[-1,0]"));
    AffineFunction<R> dynamic2(Matrix<R>("[-2,-1; 1,-2]"),Vector<R>("[1,0]"));
    //AffineFunction<R> dynamic2(Matrix<R>("[0,0;0,0]"),Vector<R>("[-1,-2]"));
    PolyhedralSet<R> activationset11(Box<R>("[-7.5,7.5]x[-3,-2]"));
    PolyhedralSet<R> activationset21(Box<R>("[-7.5,7.5]x[-7,-6]"));
    PolyhedralSet<R> activationset12(Box<R>("[-7.5,7.5]x[6,7]"));
    AffineFunction<R> reset11(Matrix<R>("[1,0;0,-1]"),Vector<R>("[0,0]"));
    AffineFunction<R> reset21(Matrix<R>("[1,0;0,-1]"),Vector<R>("[0,0]"));
    AffineFunction<R> reset12(Matrix<R>("[1,0;0,-1]"),Vector<R>("[0,0]"));
    
    AffineFunction<R> space(Matrix<R>("[1,0]"),Vector<R>("[-7.5]"));
    AffineFunction<R> invariant1(Matrix<R>("[1,0]"),Vector<R>("[-7.5]"));
    AffineFunction<R> invariant2=invariant1;
    AffineFunction<R> activation11=invariant1;
    AffineFunction<R> activation21=invariant1;
    AffineFunction<R> activation12=invariant1;

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
    R maximum_enclosure_radius=0.5;
    R grid_length=0.25;
    
    EvolutionParameters<R> parameters;
    parameters.set_maximum_enclosure_radius(maximum_enclosure_radius);
    parameters.set_grid_length(grid_length);
    parameters.set_maximum_step_size(maximum_step_size);
    parameters.set_lock_to_grid_time(lock_to_grid_time);
    parameters.set_verbosity(0);
    
    StandardApplicator<ES> applicator;
    AffineIntegrator<ES> integrator;
    StandardSatisfier<ES> satisfier;
    StandardSubdivider<ES> subdivider;
    CascadeReducer<ES> reducer(3);
    hybrid_evolver=Evolver<Sys,ES> (parameters,applicator,integrator,satisfier,subdivider,reducer);
    
    HybridBasicSet< Zonotope<R> > initial_set(mode1_id,Zonotope<R>(Box<R>("[-6.96875,-6.9375]x[-6.96875,-6.9375]")));
    cout << "initial_set=" << initial_set << endl;
    
  }

  void test_lower_semantics() {  
    cout << "Computing timed evolved set" << endl;
    ListSet< HybridBasicSet< Zonotope<R> > > lower_evolve=hybrid_evolver.reach(automaton,initial_set,Rational(1),lower_semantics);
    ListSet< HybridBasicSet< Zonotope<R> > > lower_reach=hybrid_evolver.reach(automaton,initial_set,Rational(1),lower_semantics);
    HybridListSet< Zonotope<R> > lower_evolve_sets(automaton.locations(),lower_evolve);
    HybridListSet< Zonotope<R> > lower_reach_sets(automaton.locations(),lower_reach);
    cout << "Reached (" << lower_reach_sets[mode1_id].size() << "," << lower_reach_sets[mode2_id].size() << ") enclosures, "
         << "with (" << lower_evolve_sets[mode1_id].size() << "," << lower_evolve_sets[mode2_id].size() << ") at final time."
         << endl << endl;

    epsfstream eps;
    eps.open("test_hybrid_evolver-lower.eps",bounding_box.neighbourhood(0.5));
    eps << fill_colour(white) << bounding_box;
    eps << line_style(true);
    eps << fill_colour(red) << lower_reach_sets[mode1_id];
    eps << fill_colour(magenta) << lower_reach_sets[mode2_id];
    eps << fill_colour(green) << lower_evolve_sets[mode1_id];
    eps << fill_colour(cyan) << lower_evolve_sets[mode2_id];
    eps << fill_colour(blue) << initial_set.set();
    eps.close();
  }

  void test_upper_semantics() {  
    cout << "Computing timed evolved set" << endl;
    ListSet< HybridBasicSet< Zonotope<R> > > upper_evolve=hybrid_evolver.reach(automaton,initial_set,Rational(1),upper_semantics);
    ListSet< HybridBasicSet< Zonotope<R> > > upper_reach=hybrid_evolver.reach(automaton,initial_set,Rational(1),upper_semantics);
    HybridListSet< Zonotope<R> > upper_evolve_sets(automaton.locations(),upper_evolve);
    HybridListSet< Zonotope<R> > upper_reach_sets(automaton.locations(),upper_reach);
    cout << "Reached (" << upper_reach_sets[mode1_id].size() << "," << upper_reach_sets[mode2_id].size() << ") enclosures, "
         << "with (" << upper_evolve_sets[mode1_id].size() << "," << upper_evolve_sets[mode2_id].size() << ") at final time."
         << endl << endl;

    epsfstream eps;
    eps.open("test_hybrid_evolver-upper.eps",bounding_box.neighbourhood(0.5));
    eps << fill_colour(white) << bounding_box;
    eps << line_style(true);
    eps << fill_colour(red) << upper_reach_sets[mode1_id];
    eps << fill_colour(magenta) << upper_reach_sets[mode2_id];
    eps << fill_colour(green) << upper_evolve_sets[mode1_id];
    eps << fill_colour(cyan) << upper_evolve_sets[mode2_id];
    eps << fill_colour(blue) << initial_set.set();
    eps.close(); 
  }

  
  void test() {
    ARIADNE_TEST_CALL(test_lower_semantics());
    ARIADNE_TEST_CALL(test_upper_semantics());
  }

};

*/


template<class R>
class TestConstraintBasedHybridEvolver 
{  
  typedef HybridAutomaton<R> Sys;
  typedef TaylorSet<R> EnclosureSetType;
  typedef TaylorSet<R> ES;
  typedef HybridBasicSet<ES> HES;
  typedef ListSet<HES> HLS;

  HybridEvolver< R > hybrid_evolver;
  HybridAutomaton<R> automaton;
  HybridBasicSet<ES> initial_set;
  Box<R> bounding_box;
  DiscreteState mode1_id, mode2_id;

 public:
  TestConstraintBasedHybridEvolver() 
    : automaton(""),
      initial_set(DiscreteState(2),EnclosureSetType(2)),
      bounding_box(Box<R>("[-4,4]x[-4,4]")),
      mode1_id(2),
      mode2_id(3)
  {
    hybrid_evolver_verbosity = 2;
    integrator_verbosity = 0;
    
    
  
    Box<R> r("[-1,1]x[-1,1]");
    cout << "r=" << r << endl;
    
    AffineFunction<R> dynamic(Matrix<R>("[0.25,-1.00;1.00,0.25]"),Vector<R>("[0.00,0.00]"));
    cout << "dynamic=" << dynamic << endl;
    AffineFunction<R> reset(Matrix<R>("[-0.125,0;0,-0.125]"),Vector<R>("[0,0]"));
    cout << "reset=" << reset << endl;
    
    AffineFunction<R> invariant1(Matrix<R>("[-1,0]"),Vector<R>("[1]"));
    AffineFunction<R> invariant2(Matrix<R>("[-1,0]"),Vector<R>("[4]"));
    cout << "invariant1=" << invariant1 << endl;
    cout << "invariant2=" << invariant2 << endl;
    AffineFunction<R> activation12(Matrix<R>("[-1,0]"),Vector<R>("[1.125]"));
    AffineFunction<R> guard21(Matrix<R>("[1,0]"),Vector<R>("[-1]"));
    cout << "activation12=" << activation12 << endl;
    cout << "activation21=" << guard21 << endl;
    cout << endl;
    
    automaton =HybridAutomaton<R>("Constraint-based affine test automaton");
    DiscreteState dstate1(0);
    DiscreteState dstate2(1);
    const DiscreteMode<R>& mode1=automaton.new_mode(dstate1,dynamic,invariant1);
    const DiscreteMode<R>& mode2=automaton.new_mode(dstate2,dynamic,invariant2);
    DiscreteEvent event(5);
    const DiscreteTransition<R>& transition12=automaton.new_transition(event,dstate1,dstate2,reset,activation12);
    const DiscreteTransition<R>& transition21=automaton.new_transition(event,dstate2,dstate1,reset,guard21);
    
    cout << mode1 << " " << mode2 << endl;
    cout << transition21 << " " << transition12 << " " << endl;
    cout << endl;
    
    time_type maximum_step_size=0.125;
    time_type lock_to_grid_time=0.25;
    R maximum_enclosure_radius=0.5;
    R grid_length=0.25;
    
    EvolutionParameters<R> parameters;
    parameters.set_maximum_enclosure_radius(maximum_enclosure_radius);
    parameters.set_grid_length(grid_length);
    parameters.set_maximum_step_size(maximum_step_size);
    parameters.set_lock_to_grid_time(lock_to_grid_time);
    parameters.set_verbosity(0);
    
    hybrid_evolver=HybridEvolver<R> ();
    
    HybridBasicSet<EnclosureSetType> initial_set(mode1_id,EnclosureSetType(Box<R>("[-6.96875,-6.9375]x[-6.96875,-6.9375]")));
    cout << "initial_set=" << initial_set << endl;
    
  }

  void test_lower_semantics() {  
    cout << "Computing timed evolved set" << endl;
    ListSet< HybridBasicSet< EnclosureSetType > > lower_evolve=hybrid_evolver.reach(automaton,initial_set,Rational(1),lower_semantics);
    ListSet< HybridBasicSet< EnclosureSetType > > lower_reach=hybrid_evolver.reach(automaton,initial_set,Rational(1),lower_semantics);
    HybridListSet< EnclosureSetType > lower_evolve_sets(automaton.locations(),lower_evolve);
    HybridListSet< EnclosureSetType > lower_reach_sets(automaton.locations(),lower_reach);
    cout << "Reached (" << lower_reach_sets[mode1_id].size() << "," << lower_reach_sets[mode2_id].size() << ") enclosures, "
         << "with (" << lower_evolve_sets[mode1_id].size() << "," << lower_evolve_sets[mode2_id].size() << ") at final time."
         << endl << endl;

    epsfstream eps;
    eps.open("test_hybrid_evolver-lower.eps",bounding_box.neighbourhood(0.5));
    eps << fill_colour(white) << bounding_box;
    eps << line_style(true);
    eps << fill_colour(red) << lower_reach_sets[mode1_id];
    eps << fill_colour(magenta) << lower_reach_sets[mode2_id];
    eps << fill_colour(green) << lower_evolve_sets[mode1_id];
    eps << fill_colour(cyan) << lower_evolve_sets[mode2_id];
    eps << fill_colour(blue) << initial_set.set();
    eps.close();
  }

  void test_upper_semantics() {  
    cout << "Computing timed evolved set" << endl;
    ListSet< HybridBasicSet< EnclosureSetType > > upper_evolve=hybrid_evolver.reach(automaton,initial_set,Rational(1),upper_semantics);
    ListSet< HybridBasicSet< EnclosureSetType > > upper_reach=hybrid_evolver.reach(automaton,initial_set,Rational(1),upper_semantics);
    HybridListSet< EnclosureSetType > upper_evolve_sets(automaton.locations(),upper_evolve);
    HybridListSet< EnclosureSetType > upper_reach_sets(automaton.locations(),upper_reach);
    cout << "Reached (" << upper_reach_sets[mode1_id].size() << "," << upper_reach_sets[mode2_id].size() << ") enclosures, "
         << "with (" << upper_evolve_sets[mode1_id].size() << "," << upper_evolve_sets[mode2_id].size() << ") at final time."
         << endl << endl;

    epsfstream eps;
    eps.open("test_hybrid_evolver-upper.eps",bounding_box.neighbourhood(0.5));
    eps << fill_colour(white) << bounding_box;
    eps << line_style(true);
    eps << fill_colour(red) << upper_reach_sets[mode1_id];
    eps << fill_colour(magenta) << upper_reach_sets[mode2_id];
    eps << fill_colour(green) << upper_evolve_sets[mode1_id];
    eps << fill_colour(cyan) << upper_evolve_sets[mode2_id];
    eps << fill_colour(blue) << initial_set.set();
    eps.close(); 
  }

  
  void test() {
    //ARIADNE_TEST_CALL(test_lower_semantics());
    ARIADNE_TEST_CALL(test_upper_semantics());
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
  TestConstraintBasedHybridEvolver<Flt>().test();
  cerr << "INCOMPLETE ";
  return ARIADNE_TEST_FAILURES;
}

