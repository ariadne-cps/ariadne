/***************************************************************************
 *            lorenz_attractor.cc
 *
 *  Copyright  2008  Pieter Collins
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

#include "numeric/float.h"
#include "numeric/approximate_float.h"
#include "geometry/point.h"
#include "geometry/box.h"
#include "geometry/grid.h"
#include "geometry/grid_set.h"
#include "geometry/partition_tree_set.h"
#include "geometry/rectangular_set.h"
#include "geometry/zonotope.h"
#include "evaluation/evolution_parameters.h"
#include "evaluation/vector_field_evolver.h"
#include "evaluation/standard_integrator.h"
#include "evaluation/standard_subdivider.h"
#include "evaluation/orthogonal_reducer.h"
#include "evaluation/standard_approximator.h"
#include "evaluation/fast_approximator.h"
#include "evaluation/reachability_analyser.h"
#include "output/epsstream.h"
#include "output/logging.h"
#include "models/lorenz.h"

using namespace Ariadne;
using namespace Ariadne::Models;
using namespace std;

template<class R> int lorenz_attractor();

int main() {
  return lorenz_attractor<Float64>();
}

template<class R> 
int 
lorenz_attractor()
{
  typedef Zonotope<R> ES;

  // Specify the time for finite-time evolution
  Rational time(2);

  // Specify the evolution parameters
  EvolutionParameters<R> evolution_parameters;
  double maximum_enclosure_radius=0.25;
  double grid_length=0.125;
  int subdivisions=128;
  evolution_parameters.set_maximum_enclosure_radius(maximum_enclosure_radius);
  evolution_parameters.set_grid_length(grid_length);
  evolution_parameters.set_lock_to_grid_time(0.25);
  evolution_parameters.set_maximum_step_size(0.25);
  //evolution_parameters.set_bounding_domain_size(1.0);
  
  StandardIntegrator<ES> integrator;
  StandardSubdivider<ES> subdivider;
  OrthogonalReducer<ES> reducer;
  Evolver<VectorField<R>,ES> evolver(evolution_parameters,integrator,subdivider,reducer);

  //StandardApproximator<ES> approximator;
  FastApproximator<ES> approximator;
  
  set_applicator_verbosity(0);
  
  ReachabilityAnalyser< VectorField<R>, GridApproximationScheme<R> > analyser(evolution_parameters,evolver,approximator);

  double parameter_array[]={8./3.,10,28};
  Point<R> parameters=Point<R>(3,parameter_array);
  R beta=parameters[0];
  R rho=parameters[1];
  R sigma=parameters[2];


  LorenzSystem<R> lorenz=LorenzSystem<R>(parameters);
  cout << "system = " << lorenz << endl;
  Point<R> pt("(1.0, 0.0, 0.0)"); 
  smoothness_type s=3;
  cout << "pt="<<pt<<" s="<<s<<endl;
  cout << "lorenz(pt)="<<lorenz(pt)<<endl;
  cout << "lorenz.jacobian(pt)="<<lorenz.jacobian(pt) << endl; 
  cout << "lorenz.derivative(pt,s)="<<lorenz.derivative(pt,s) << endl; 

  Box<R> bounding_box("[-8.0,8.0]x[-8.0,8.0]x[-8.0,8.0]") ;
  Box<R> initial_box("[0.9,1.1]x[-0.1,0.1]x[-0.1,0.1]"); // initial state
  Box<R> graphic_bounding_box=Box<R>("[-8.1,8.1]x[-8.1,8.1]x[-8.1,8.1]"); // eps bounding box
  
  Zonotope<R> initial_zonotope(initial_box);
  ImageSet<R> initial_set(initial_box);

  cout << "Computing finite time evolution..." << flush;
  ListSet< Zonotope<R> > evolve_set=evolver.evolve(lorenz,initial_zonotope,time);
  cout << "  done" << endl;
  cout << "Computing finite time reachable set..." << flush;
  ListSet< Zonotope<R> > reach_set=evolver.reach(lorenz,initial_zonotope,time);
  cout << "  done" << endl;

  
  cout << "Computing attractor..." << flush;
  SetInterface< Box<R> >* cr=analyser.chain_reach(lorenz,initial_set);
  cout << "  done" << endl;
  GridMaskSet<R> grid_chain_reach_set=*dynamic_cast<GridMaskSet<R>*>(cr);
  PartitionTreeSet<R> tree_chain_reach_set(grid_chain_reach_set);
   

  /*
  cout << "Skipping computation of attractor..." << flush;
  GridMaskSet<R> grid_chain_reach_set(Grid<R>(3,0.125),bounding_box);
  PartitionTreeSet<R> tree_chain_reach_set(grid_chain_reach_set);
  */

  cout << grid_chain_reach_set <<endl;
  cout << tree_chain_reach_set <<endl;
  
  cout << "evolve_set.size()=" << evolve_set.size() << endl;
  cout << "reach_set.size()=" << reach_set.size() << endl;
  cout << "grid_chain_reach_set.size()=" << grid_chain_reach_set.size() << endl;
  cout << "tree_chain_reach_set.size()=" << tree_chain_reach_set.size() << "  " << flush;
  cout << "tree_chain_reach_set.capacity()=" << tree_chain_reach_set.capacity() << endl;

  epsfstream eps;
  eps.open("lorenz-reach.eps",graphic_bounding_box);
  eps << fill_colour(white) << bounding_box;
  eps << fill_colour(green) << reach_set;
  eps << fill_colour(yellow) << evolve_set;
  eps << fill_colour(blue) << initial_box;
  eps.close();

  eps.open("lorenz-chainreach.eps",graphic_bounding_box);
  eps << fill_colour(white) << bounding_box;
  eps << line_style(false);
  eps << fill_colour(green) << tree_chain_reach_set;
  eps << fill_colour(blue) << initial_box;
  eps << line_style(true);
  eps << tree_chain_reach_set.partition_tree();
  eps.close();


  return 0;
}
