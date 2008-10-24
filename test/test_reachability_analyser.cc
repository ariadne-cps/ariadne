/***************************************************************************
 *            test_reachability_analyser.cc
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

#include "approximate_taylor_model.h"
#include "grid_set.h"
#include "hybrid_set.h"
#include "hybrid_automaton.h"
#include "evolution_parameters.h"
#include "hybrid_evolver.h"
#include "discretiser.h"
#include "analyser.h"
#include "graphics.h"
#include "logging.h"

#include "models.h"

#include "test.h"

using namespace std;
using namespace Ariadne;
using Ariadne::Models::Henon;



class TestReachabilityAnalyser 
{  

  HybridAnalyser analyser;
  HybridAutomaton system;
  Grid grid;
  Interval bound;
  HybridImageSet initial_set;
  HybridTime reach_time;
 
  typedef ApproximateTaylorModel EnclosureType;

 public:
  static HybridAnalyser build_analyser()
  {
    EvolutionParameters parameters;
    parameters.maximum_enclosure_radius=0.25;
    parameters.grid_length=0.0625;
    parameters.maximum_step_size=0.125;
    parameters.lock_to_grid_time=0.5;
    
    Grid grid(2);
    HybridEvolver evolver(parameters);
    EvolverInterface<HybridAutomaton,DefaultHybridEnclosureType>& evolver_interface
      =evolver;
    //HybridDiscretiser<EnclosureType> discretiser(evolver);
    return HybridAnalyser(parameters,evolver_interface);
  }

  TestReachabilityAnalyser()
    : analyser(build_analyser()),
      system("Henon Map"),
      grid(make_vector<Float>("[0.25,0.25]")),
      bound(-8,8),
      reach_time(2.0,3)
  {
    DiscreteState loc(1);

    Function<Henon> henon(make_point("(1.5,-0.375)"));
    system.new_mode(DiscreteState(1),ConstantFunction(Vector<Float>(0),2));
    system.new_forced_transition(DiscreteEvent(1),DiscreteState(1),DiscreteState(1),henon,ConstantFunction(Vector<Float>(1,1.0),2));

    ImageSet initial_box(make_box("[-6.96875,-6.9375]x[-6.96875,-6.9375]"));
    initial_set[loc]=initial_box;
    cout << "initial_set=" << initial_set << endl;

    //ARIADNE_ASSERT(subset(initial_set[loc],bounding_set[loc]));
    
  }

  template<class S, class IS> void plot(const char* name, const Box& bounding_box, const S& set, const IS& initial_set) {
    Graphic g;
    g << fill_colour(white) << bounding_box;
    g << line_style(true);
    g << fill_colour(red) << set;
    g << fill_colour(blue);
    g << initial_set;
    g.write(name);
  }

  void test_lower_reach() {  
    cout << "test_lower_reach() not used\n";
  /*
    cout << "Computing timed reachable set" << endl;
    GridMaskSet lower_reach=map_analyser->lower_reach(system,initial_set,Rational(1));
    cout << "Reached " << lower_reach.size() << " cells out of " << lower_reach.capacity() << endl << endl;
    plot("test_reachability_analyser-map_lower_reach.eps",bounding_box,lower_reach,initial_set);
  */
  }
  
  void test_upper_reach() {  
    cout << "Computing timed reachable set" << endl;
    DiscreteState loc(1);
    Box bounding_box(2,bound);
    HybridGridTreeSet& upper_reach_set=*analyser.upper_reach(system,initial_set,reach_time);
    const GridTreeSet& upper_reach=upper_reach_set[loc];
    ImageSet& initial=initial_set[loc];
    //cout << "Reached " << upper_reach.size() << " cells out of " << upper_reach.capacity() << endl << endl;
    plot("test_reachability_analyser-map_upper_reach.eps",bounding_box,upper_reach,initial);
  }
  
  void test_chain_reach() {  
    cout << "Computing chain reachable set" << endl;
    DiscreteState loc(1);
    HybridBoxes bounding_boxes
      =Ariadne::bounding_boxes(system.state_space(),bound);
    Box bounding_box=bounding_boxes[loc];

    HybridGridTreeSet& chain_reach_set=*analyser.chain_reach(system,initial_set,bounding_boxes);
    plot("test_reachability_analyser-map_chain_reach.eps",bounding_box,chain_reach_set[loc],initial_set[loc]);
  }
  
  void test() {
    ARIADNE_TEST_CALL(test_lower_reach());
    ARIADNE_TEST_CALL(test_upper_reach());
    ARIADNE_TEST_CALL(test_chain_reach());
  }

};


int main(int nargs, const char* args[]) 
{
  TestReachabilityAnalyser().test();
  cerr << "INCOMPLETE ";
  //cerr << "SKIPPED ";
  //++ARIADNE_TEST_FAILURES;
  return ARIADNE_TEST_FAILURES;
}

