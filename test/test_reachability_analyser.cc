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

#include "ariadne.h"
#include "test_float.h"
#include "geometry/set_interface.h"
#include "geometry/set_reference.h"
#include "geometry/hybrid_set.h"
#include "geometry/zonotope.h"
#include "geometry/empty_set.h"
#include "geometry/polyhedral_set.h"
#include "system/transition_system.h"
#include "system/numerical_system.h"
#include "system/affine_map.h"
#include "system/affine_vector_field.h"
#include "system/hybrid_automaton.h"
#include "evaluation/evolution_parameters.h"
#include "evaluation/standard_applicator.h"
#include "evaluation/standard_subdivider.h"
#include "evaluation/cascade_reducer.h"
#include "evaluation/map_evolver.h"
#include "evaluation/standard_approximator.h"
#include "evaluation/discretiser.h"
#include "evaluation/reachability_analyser.h"
#include "output/epsstream.h"
#include "output/logging.h"
#include "models/henon.h"

#include "test.h"

using namespace std;
using namespace Ariadne;
using Ariadne::Models::HenonMap;


template<class R>
class TestReachabilityAnalyser 
{  
  typedef Map<R> Sys;
  typedef GridApproximationScheme<R> Aprx;
  typedef Zonotope<R> ES;

  ReachabilityAnalyser<Sys,Aprx> map_analyser;
  Map<R> system;
  Grid<R> grid;
  Box<R> bounding_box;
  GridMaskSet<R> initial_set;
  GridMaskSet<R> bounding_set;

 public:
  static ReachabilityAnalyser< Map<R>, GridApproximationScheme<R> > build_map_analyser()
  {
    EvolutionParameters<R> parameters;
    parameters.set_maximum_basic_set_radius(0.25);
    parameters.set_grid_length(0.0625);
    parameters.set_maximum_step_size(0.125);
    parameters.set_lock_to_grid_time(0.5);
    parameters.set_verbosity(0);
    
    Grid<R> grid(2,parameters.grid_length());
    StandardApplicator<ES> applicator;
    StandardSubdivider<ES> subdivider;
    CascadeReducer<ES> reducer(3);
    Evolver<Sys,ES> map_evolver(parameters,applicator,subdivider,reducer);
    StandardApproximator<ES> approximator(grid);
    Discretiser<Sys,Aprx,ES> discretiser(parameters,map_evolver,approximator);
    return ReachabilityAnalyser<Sys,Aprx>(parameters,discretiser);
  }

  TestReachabilityAnalyser()
    : map_analyser(build_map_analyser()),
      system(HenonMap<R>(Point<R>("(1.5,-0.375)"))),
      grid(Vector<R>("[0.25,0.25]")),
      bounding_box("[-8,8]x[-8,8]"),
      initial_set(grid,bounding_box),
      bounding_set(grid,bounding_box)
  {
    RectangularSet<R> initial_rectangle(Box<R>("[-6.96875,-6.9375]x[-6.96875,-6.9375]"));

    initial_set.adjoin_outer_approximation(initial_rectangle);
    cout << "initial_set=" << this->initial_set << endl;
    bounding_set.adjoin_over_approximation(bounding_box);

    ARIADNE_ASSERT(subset(initial_set,bounding_set));
    
  }

  template<class S, class IS> void plot(const char* name, const Box<R>& bounding_box, const S& set, const IS& intial_set) {
    epsfstream eps;
    eps.open(name,bounding_box.neighbourhood(0.5));
    eps << fill_colour(white) << bounding_box;
    eps << line_style(true);
    eps << fill_colour(red) << set;
    eps << fill_colour(blue) << initial_set;
    eps.close();
  }

  void test_lower_reach() {  
    cout << "test_lower_reach() not used\n";
  /*
    cout << "Computing timed reachable set" << endl;
    GridMaskSet<R> lower_reach=map_analyser->lower_reach(system,initial_set,Rational(1));
    cout << "Reached " << lower_reach.size() << " cells out of " << lower_reach.capacity() << endl << endl;
    plot("test_reachability_analyser-map_lower_reach.eps",bounding_box,lower_reach,initial_set);
  */
  }
  
  void test_upper_reach() {  
    cout << "Computing timed reachable set" << endl;
    SetInterface< Box<R> >* upper_reach_set=map_analyser.upper_reach(system,initial_set,Integer(1));
    GridMaskSet<R> upper_reach=*dynamic_cast<GridMaskSet<R>*>(upper_reach_set);
    cout << "Reached " << upper_reach.size() << " cells out of " << upper_reach.capacity() << endl << endl;
    plot("test_reachability_analyser-map_upper_reach.eps",bounding_box,upper_reach,initial_set);
  }
  
  void test_chain_reach() {  
    cout << "Computing chain reachable set" << endl;
    SetInterface< Box<R> >* chain_reach_set=map_analyser.chain_reach(system,initial_set);
    GridMaskSet<R> chain_reach=*dynamic_cast<GridMaskSet<R>*>(chain_reach_set);
    cout << "Reached " << chain_reach.size() << " cells out of " << chain_reach.capacity() << endl << endl;
    plot("test_reachability_analyser-map_chain_reach.eps",bounding_box,chain_reach,initial_set);
  }
  
  void test() {
    ARIADNE_TEST_CALL(test_lower_reach());
    ARIADNE_TEST_CALL(test_upper_reach());
    ARIADNE_TEST_CALL(test_chain_reach());
  }

};


int main(int nargs, const char* args[]) 
{
  TestReachabilityAnalyser<Flt>().test();
  cerr << "INCOMPLETE ";
  //cerr << "SKIPPED ";
  //++ARIADNE_TEST_FAILURES;
  return ARIADNE_TEST_FAILURES;
}

