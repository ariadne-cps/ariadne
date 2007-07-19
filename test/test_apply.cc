/***************************************************************************
 *            test_apply.cc
 *
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
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
#include "evaluation/applicator.h"
#include "output/epsfstream.h"
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

template<class R> int test_apply();

int main() {
  return test_apply<Float>();
}

template<class R> 
int 
test_apply()
{
  //set_evaluation_verbosity(0);
  typedef Interval<R> I;

  double basic_set_radius=0.25;
  double grid_size=0.5;
  double bound=4.0;

  R a=1.5;
  R b=0.5;

  HenonMap<R> henon=HenonMap<R>(a,b);
  HenonInverseMap<R> henon_inverse=HenonInverseMap<R>(a,b);
  cout << henon << endl << henon_inverse << endl;

  Rectangle<R> bounding_box=Rectangle<R>("[-4,4]x[-4,4]") ;
  Rectangle<R> eps_bounding_box=bounding_box.neighbourhood(0.1);
  
  Grid<R> grid(2,grid_size);
  FiniteGrid<R> finite_grid=FiniteGrid<R>(grid,bounding_box); // grid

  Rectangle<R> r=Rectangle<R>("[1.499,1.501]x[0.499,0.501]"); // initial state
  Zonotope<R,R> z=Zonotope<R,R>(Rectangle<R>("[1.499,1.501]x[0.499,0.501]")); // initial state
  Zonotope<I,R> ez=Zonotope<I,R>(Rectangle<R>("[1.499,1.501]x[0.499,0.501]")); // initial state
  Zonotope<I,I> iz=Zonotope<I,I>(Rectangle<R>("[1.499,1.501]x[0.499,0.501]")); // initial state
  Polytope<R> pl=Polytope<R>(Rectangle<R>("[1.499,1.501]x[0.499,0.501]")); // initial state
  
  Applicator<R> apply(basic_set_radius,grid_size);

  apply.set_grid_size(grid_size);
  apply.set_maximum_basic_set_radius(4.0*grid_size);
  apply.set_default_bound(bound);
  

  //Test evaluation on different classes of sets
  Rectangle<R> fr=apply.evaluate(henon,r);
  Zonotope<R,R> fz=apply.evaluate(henon,z);
  Zonotope<I,R> fez=apply.evaluate(henon,ez);
  Zonotope<I,I> fiz=apply.evaluate(henon,iz);

  Rectangle<R> pfr=apply.evaluate(henon_inverse,fr);
  cout << "r=" << r << " fr=" << fr << " pfr="<< pfr << endl;
  Zonotope<I,R> pfez=apply.evaluate(henon_inverse,fez);
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

  GridMaskSet<R> grid_image_set=apply.image(henon,grid_initial_set,grid_bounding_set);
  cout << "grid_image_set=" << grid_image_set << endl;
  GridMaskSet<R> grid_preimage_set=apply.preimage(henon,grid_initial_set,grid_bounding_set);
  cout << "grid_preimage_set=" << grid_preimage_set << endl;
  GridMaskSet<R> grid_inverse_image_set=apply.image(henon_inverse,grid_initial_set,grid_bounding_set);
  cout << "grid_inverse_image_set=" << grid_inverse_image_set << endl;

  ListSet< Zonotope<I,R> > list_initial_set(grid_initial_set);
  cout << "list_initial_set=" << list_initial_set << endl;
  ListSet< Zonotope<I,R> > list_reach_set=apply.reach(henon,list_initial_set);
  cout << "list_reach_set=" << list_reach_set << endl;

  GridMultiMap<R> discretization=apply.discretize(henon,grid_bounding_set,grid);
  cout << "discretization=GridMultiMap(...)" << endl;

  epsfstream eps;
  eps.open("test_apply-1.eps",eps_bounding_box);
  eps.set_fill_colour("green");
  eps<<list_reach_set;
  eps.set_fill_colour("blue");
  eps<<grid_initial_set;
  eps.set_fill_colour("green");
  eps<<grid_image_set;
  eps.set_fill_colour("red");
  eps<<grid_inverse_image_set;
  eps.set_fill_colour("yellow");
  eps<<grid_preimage_set;
  eps.set_fill_colour("magenta");
  eps<<list_reach_set;
  eps.close();

  eps.open("test_apply-2.eps",eps_bounding_box);
  eps.set_line_style(true);
  eps.set_fill_colour("green");
  eps<<grid_image_set;
  eps.set_fill_colour("red");
  eps<<apply.preimage(henon,grid_image_set,grid_bounding_set);
  eps.set_fill_colour("blue");
  eps<<grid_initial_set;
  eps.close();


  cout << endl;

  cout << "Computing with SetInterface" << endl;

  cout << "Computing image set" << endl;
  shared_ptr< SetInterface<R> > image_set_ptr(apply.image(henon,initial_set));
  cout << "image_set=" << *image_set_ptr << endl;
  cout << "Computing preimage set" << endl;
  apply.set_grid_size(grid_size/2);
  shared_ptr< SetInterface<R> > preimage_set_ptr(apply.preimage(henon,initial_set,bounding_set));
  apply.set_grid_size(grid_size);
  cout << "preimage_set=" << *preimage_set_ptr << endl;
  cout << "Computing inverse image set" << endl;
  shared_ptr< SetInterface<R> > inverse_image_set_ptr(apply.image(henon_inverse,initial_set));
  cout << "inverse_image_set=" << *inverse_image_set_ptr << endl;

  eps.open("test_apply-3.eps",eps_bounding_box);
  eps.set_fill_colour("blue");
  eps<<initial_set;
  eps.set_fill_colour("green");
  eps<<*image_set_ptr;
  eps.set_fill_colour("red");
  eps<<*inverse_image_set_ptr;
  eps.set_fill_colour("yellow");
  eps<<*preimage_set_ptr;
  eps.close();
  

  cout << "Computing reach set" << endl;
  apply.set_grid_size(0.25*grid_size);
  shared_ptr< SetInterface<R> > reach_set_ptr(apply.reach(henon,initial_set));
  cout << "reach_set=" << *reach_set_ptr << endl;
  apply.set_grid_size(grid_size);
  cout << "Computing chainreach set" << endl;
  shared_ptr< SetInterface<R> > chainreach_set_ptr(apply.chainreach(henon,initial_set,bounding_set));
  cout << "chainreach_set=" << *chainreach_set_ptr << endl;

  eps.open("test_apply-4.eps",eps_bounding_box);
  eps.set_fill_colour("cyan");
  eps << bounding_set;
  eps.set_line_style(false);
  eps.set_fill_colour("green");
  eps << *chainreach_set_ptr;
  eps.set_line_style(true);
  eps.set_fill_colour("magenta");
  eps << *reach_set_ptr;
  eps.set_fill_style(false);
  eps << initial_set;
  eps << bounding_set;
  eps.close();
  

  cout << "Computing viability kernel" << endl;
  apply.set_grid_size(grid_size);
  shared_ptr< SetInterface<R> > viability_kernel_ptr(apply.viable(henon,bounding_set));
  cout << "viability_kernel=" << *viability_kernel_ptr << endl;

  eps.open("test_apply-5.eps",eps_bounding_box);
  eps.set_fill_style(true);
  eps.set_fill_colour("cyan");
  eps << bounding_set;
  eps.set_fill_colour("magenta");
  eps << *viability_kernel_ptr;
  eps.close();

  return 0;
}
