/***************************************************************************
 *            test_map_evolver.cc
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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

#include "test_float.h"

#include "base/pointer.h"
#include "geometry/point.h"
#include "geometry/zonotope.h"
#include "geometry/polytope.h"
#include "geometry/rectangular_set.h"
#include "system/grid_multimap.h"
#include "evaluation/evolution_parameters.h"
#include "evaluation/map_evolver.h"
#include "evaluation/standard_applicator.h"
#include "output/epsstream.h"
#include "output/logging.h"

#include "models/henon.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::Combinatoric;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;
using namespace Ariadne::Output;
using namespace std;

template<class R> int test_map_evolver();

int main() {
  return test_map_evolver<Flt>();
}

template<class R> 
int 
test_map_evolver()
{
  set_evaluation_verbosity(0);
  typedef Interval<R> I;
  typedef Zonotope<R> ZBS;
  typedef Polytope<R> PBS;

  Integer maximum_number_of_steps=100;
  R maximum_basic_set_radius=0.25;
  R grid_length=0.5;
  R fine_grid_length=0.5/16;
  R bounding_domain_size=4.0;

  R a=1.5;
  R b=0.5;
  R p[2]={a,b};
  
  HenonMap<R> henon=HenonMap<R>(Point<R>(2,p));
  HenonInverseMap<R> henon_inverse=HenonInverseMap<R>(henon.parameters());
  cout << henon << endl << henon_inverse << endl;

  Box<R> bounding_box=Box<R>("[-4,4]x[-4,4]") ;
  Box<R> eps_bounding_box=bounding_box.neighbourhood(0.1);
  
  Grid<R> grid(2,grid_length);
  Grid<R> fine_grid(2,fine_grid_length);
  FiniteGrid<R> finite_grid=FiniteGrid<R>(grid,bounding_box); // grid

  Box<R> bx=Box<R>("[1.499,1.501]x[0.499,0.501]");
  Zonotope<R> z(bx);
  Polytope<R> pl(bx);
  
  StandardApplicator<R> standard;

  //Test evaluation on different classes of sets
  Zonotope<R> fz=standard.apply(henon,z);

  EvolutionParameters<R> parameters;
  parameters.set_maximum_basic_set_radius(maximum_basic_set_radius);
  parameters.set_grid_length(grid_length);
  parameters.set_bounding_domain_size(bounding_domain_size);

  MapEvolver<ZBS> evolver(parameters);
  Zonotope<R> pfz=standard.apply(henon_inverse,fz);
  cout << "z=" << z << " fz=" << fz << " pfz="<< pfz << endl;
  
  ConstraintSet<R> bounding_set(bounding_box);
  Zonotope<R> initial_zonotope(Box<R>("[-0.2,0.8]x[0.3,1.7]"));
  ConstraintSet<R> initial_set(Box<R>("[-0.2,0.8]x[0.3,1.7]"));
  

  // Test evaluation on concrete sets
  GridMaskSet<R> grid_initial_set(grid,bounding_box);
  grid_initial_set.adjoin_outer_approximation(initial_set);
  cout << "grid_initial_set=" << grid_initial_set << endl;
  GridMaskSet<R> grid_bounding_set(grid,bounding_box);
  grid_bounding_set.adjoin_outer_approximation(bounding_set);
  cout << "grid_bounding_set=" << grid_bounding_set << endl;

  epsfstream eps;

  cout << "Computing with SetInterface" << endl;
  
  cout << "Computing reach set" << endl;
  evolver.parameters().set_grid_length(div_approx(grid_length,4));
  shared_ptr< SetInterface<R> > reach_set_ptr(evolver.lower_reach(henon,initial_set,maximum_number_of_steps));
  cout << "reach_set=" << *reach_set_ptr << endl;
  evolver.parameters().set_grid_length(grid_length);
  cout << "Computing chainreach set" << endl;
  shared_ptr< SetInterface<R> > chainreach_set_ptr(evolver.chainreach(henon,initial_set));
  cout << "chainreach_set=" << *chainreach_set_ptr << endl;

  eps.open("test_map_evolver-abstract_reach.eps",eps_bounding_box);
  eps << fill_colour(cyan) << bounding_set;
  eps << line_style(false) << fill_colour(green) << *chainreach_set_ptr;
  eps << line_style(true) << fill_colour(magenta) << *reach_set_ptr;
  eps << fill_colour(transparant);
  eps << initial_set;
  eps << bounding_set;
  eps.close();
  
 /*
  cout << "Computing viability kernel" << endl;
  evolver.parameters().set_grid_length(grid_length);
  shared_ptr< SetInterface<R> > viability_kernel_ptr(evolver.viable(henon,bounding_set));
  cout << "viability_kernel=" << *viability_kernel_ptr << endl;

  eps.open("test_map_evolver-abstract_viable.eps",eps_bounding_box);
  eps << fill_colour(cyan) << bounding_set;
  eps << fill_colour(magenta) << *viability_kernel_ptr;
  eps.close();
 */

  return ARIADNE_TEST_FAILURES;
}
