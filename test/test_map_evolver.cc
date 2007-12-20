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
#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"
#include "geometry/polytope.h"
#include "geometry/rectangular_set.h"
#include "system/grid_multimap.h"
#include "evaluation/evolution_parameters.h"
#include "evaluation/map_evolver.h"
#include "evaluation/applicator.h"
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
  typedef Zonotope<I,R> BS;

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

  Rectangle<R> bounding_box=Rectangle<R>("[-4,4]x[-4,4]") ;
  Rectangle<R> eps_bounding_box=bounding_box.neighbourhood(0.1);
  
  Grid<R> grid(2,grid_length);
  Grid<R> fine_grid(2,fine_grid_length);
  FiniteGrid<R> finite_grid=FiniteGrid<R>(grid,bounding_box); // grid

  Rectangle<R> r=Rectangle<R>("[1.499,1.501]x[0.499,0.501]"); // initial state
  Zonotope<R> z=Zonotope<R>(Rectangle<R>("[1.499,1.501]x[0.499,0.501]")); // initial state
  Zonotope<R,UniformErrorTag> ez=Zonotope<R,UniformErrorTag>(Rectangle<R>("[1.499,1.501]x[0.499,0.501]")); // initial state
  Polytope<R> pl=Polytope<R>(Rectangle<R>("[1.499,1.501]x[0.499,0.501]")); // initial state
  
  //Test evaluation on different classes of sets
  Rectangle<R> fr=apply(henon,r);
  Zonotope<R> fz=apply(henon,z);
  Zonotope<R,UniformErrorTag> fez=apply(henon,ez);

  EvolutionParameters<R> parameters;
  parameters.set_maximum_basic_set_radius(maximum_basic_set_radius);
  parameters.set_grid_length(grid_length);
  parameters.set_bounding_domain_size(bounding_domain_size);

  MapEvolver<R> evolver(parameters);
  Rectangle<R> pfr=apply(henon_inverse,fr);
  cout << "r=" << r << " fr=" << fr << " pfr="<< pfr << endl;
  Zonotope<R,UniformErrorTag> pfez=apply(henon_inverse,fez);
  cout << "ez=" << ez << " fez=" << fez << " pfez="<< pfez << endl;
  
  RectangularSet<R> bounding_set(bounding_box);
  RectangularSet<R> initial_set("[-0.2,0.8]x[0.3,1.7]");
  

  // Test evaluation on concrete sets
  GridMaskSet<R> grid_initial_set(grid,initial_set);
  grid_initial_set.adjoin_over_approximation(Rectangle<R>(initial_set));
  cout << "grid_initial_set=" << grid_initial_set << endl;
  GridMaskSet<R> grid_bounding_set(grid,bounding_set);
  grid_bounding_set.adjoin_over_approximation(Rectangle<R>(bounding_set));
  cout << "grid_bounding_set=" << grid_bounding_set << endl;

  cout << endl << endl;

  cout << "Computing with concrete sets" << endl;

  GridMaskSet<R> grid_image_set=evolver.image(henon,grid_initial_set,grid_bounding_set);
  cout << "grid_image_set=" << grid_image_set << endl;
  GridMaskSet<R> grid_preimage_set=evolver.preimage(henon,grid_initial_set,grid_bounding_set);
  cout << "grid_preimage_set=" << grid_preimage_set << endl;
  GridMaskSet<R> grid_inverse_image_set=evolver.image(henon_inverse,grid_initial_set,grid_bounding_set);
  cout << "grid_inverse_image_set=" << grid_inverse_image_set << endl;
  GridMaskSet<R> grid_chainreach_set=evolver.chainreach(henon,grid_initial_set,grid_bounding_set);
  cout << "grid_chainreach_set=" << grid_chainreach_set << endl;

  RectangularSet<R> non_fixed_initial_set("[-0.5,-0.25]x[-0.625,-0.375]");
  ListSet< Rectangle<R> > list_initial_set=Geometry::point_approximation(non_fixed_initial_set,fine_grid);
  cout << "list_initial_set=" << list_initial_set << endl;
  ListSet< Rectangle<R> > list_reach_set=evolver.lower_reach(henon,list_initial_set);
  cout << "list_reach_set=" << list_reach_set << endl;
  GridMaskSet<R> list_grid_reach_set(grid,bounding_set);
  list_grid_reach_set.adjoin_outer_approximation(list_reach_set);
  cout << "list_grid_reach_set=" << list_grid_reach_set << endl;

  GridMultiMap<R> discretization=evolver.discretize(henon,grid_bounding_set,grid);
  cout << "discretization=GridMultiMap(...)" << endl;

  epsfstream eps;
  eps.open("test_map_evolver-list.eps",eps_bounding_box);
  eps << fill_colour(green) << grid_chainreach_set;
  eps << line_style(false);
  eps << fill_colour(red) << list_grid_reach_set;
  eps << fill_colour(black) << list_reach_set;
  eps << fill_colour(yellow) << list_initial_set;
  eps.close();

  eps.open("test_map_evolver-grid.eps",eps_bounding_box);
  eps << fill_colour(blue) << grid_initial_set;
  eps << fill_colour(green) << grid_image_set;
  eps << fill_colour(red) << grid_inverse_image_set;
  eps << fill_colour(yellow) << grid_preimage_set;
  eps.close();

  eps.open("test_map_evolver-grid_preimage.eps",eps_bounding_box);
  eps << line_style(true);
  eps << fill_colour(green) << grid_image_set;
  eps << fill_colour(red) << evolver.preimage(henon,grid_image_set,grid_bounding_set);
  eps << fill_colour(blue) << grid_initial_set;
  eps.close();


  cout << endl;

  cout << "Computing with SetInterface" << endl;

  cout << "Computing image set" << endl;
  shared_ptr< SetInterface<R> > image_set_ptr(evolver.image(henon,initial_set));
  cout << "image_set=" << *image_set_ptr << endl;
  cout << "Computing preimage set" << endl;
  evolver.parameters().set_grid_length(div_approx(grid_length,2));
  shared_ptr< SetInterface<R> > preimage_set_ptr(evolver.preimage(henon,initial_set,bounding_set));
  evolver.parameters().set_grid_length(grid_length);
  cout << "preimage_set=" << *preimage_set_ptr << endl;
  cout << "Computing inverse image set" << endl;
  shared_ptr< SetInterface<R> > inverse_image_set_ptr(evolver.image(henon_inverse,initial_set));
  cout << "inverse_image_set=" << *inverse_image_set_ptr << endl;

  eps.open("test_map_evolver-abstract_image.eps",eps_bounding_box);
  eps << fill_colour(blue) << initial_set;
  eps << fill_colour(green) << *image_set_ptr;
  eps << fill_colour(red) << *inverse_image_set_ptr;
  eps << fill_colour(yellow) << *preimage_set_ptr;
  eps.close();
  

  cout << "Computing reach set" << endl;
  evolver.parameters().set_grid_length(div_approx(grid_length,4));
  shared_ptr< SetInterface<R> > reach_set_ptr(evolver.lower_reach(henon,initial_set));
  cout << "reach_set=" << *reach_set_ptr << endl;
  evolver.parameters().set_grid_length(grid_length);
  cout << "Computing chainreach set" << endl;
  shared_ptr< SetInterface<R> > chainreach_set_ptr(evolver.chainreach(henon,initial_set,bounding_set));
  cout << "chainreach_set=" << *chainreach_set_ptr << endl;

  eps.open("test_map_evolver-abstract_reach.eps",eps_bounding_box);
  eps << fill_colour(cyan) << bounding_set;
  eps << line_style(false) << fill_colour(green) << *chainreach_set_ptr;
  eps << line_style(true) << fill_colour(magenta) << *reach_set_ptr;
  eps << fill_colour(transparant);
  eps << initial_set;
  eps << bounding_set;
  eps.close();
  

  cout << "Computing viability kernel" << endl;
  evolver.parameters().set_grid_length(grid_length);
  shared_ptr< SetInterface<R> > viability_kernel_ptr(evolver.viable(henon,bounding_set));
  cout << "viability_kernel=" << *viability_kernel_ptr << endl;

  eps.open("test_map_evolver-abstract_viable.eps",eps_bounding_box);
  eps << fill_colour(cyan) << bounding_set;
  eps << fill_colour(magenta) << *viability_kernel_ptr;
  eps.close();

  return 0;
}
