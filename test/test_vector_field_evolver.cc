/***************************************************************************
 *            test_integrate.cc
 *
 *  Copyright  2006  Pieter Collins
 *  Email Pieter.Collins@cwi.nl, casagrande@dimi.uniud.it
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

#include <fstream>

#include "test_float.h"

#include "ariadne.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"
#include "geometry/list_set.h"
#include "system/affine_vector_field.h"
#include "evaluation/evolution_parameters.h"
#include "evaluation/vector_field_evolver.h"
#include "evaluation/lohner_integrator.h"
#include "evaluation/affine_integrator.h"
#include "models/vanderpol.h"
#include "output/epsstream.h"
#include "output/logging.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;
using namespace Ariadne::Output;
using namespace std;

template<class R> int test_integrator();
template<class R> int test_vector_field_evolver();

int main() {
  test_integrator<Float>();
  test_vector_field_evolver<Float>();
  return 0;
}

template<class R> 
int 
test_integrator()
{
  cout << __PRETTY_FUNCTION__ << endl;
  typedef Interval<R> I;

  EvolutionParameters<R> parameters;
  parameters.set_maximum_basic_set_radius(0.125);
  parameters.set_lock_to_grid_time(0.5);
  parameters.set_maximum_step_size(0.125);
  
  // Test constructor/destructor
  LohnerIntegrator<R>* lohner_ptr;
  lohner_ptr=new LohnerIntegrator<R>();
  delete lohner_ptr;
  
  AffineVectorField<R> avf=AffineVectorField<R>(Matrix<R>("[-0.25,-1.0;+1.0,-0.25]"),Vector<R>("[0.25,0.0]"));
  cout << "avf=" << avf << endl;
  VanDerPolEquation<R> vdp=VanDerPolEquation<R>(R(0.865));
  cout << "vpd=" << avf << endl;

  Rectangle<R> r=Rectangle<R>("[0.98,1.02]x[0.48,0.52]");
  cout << "r=" << r << endl;
  Zonotope<I,R> ez=Zonotope<I,R>(r);
  cout << "ez=" << ez << endl;
  Zonotope<I,I> iz=Zonotope<I,I>(r);
  cout << "iz=" << iz << endl;

  ListSet< Zonotope<I,R> > ezls=ListSet< Zonotope<I,R> >(ez);
  ezls.adjoin(Zonotope<R>(Rectangle<R>("[1.02,1.06]x[0.48,0.52]")));
  cout << "ezls.size()=" << ezls.size() << endl;
  
  ListSet< Zonotope<I,I> > izls=ListSet< Zonotope<I> >(iz);
  izls.adjoin(Zonotope<I,I>(Rectangle<R>("[1.02,1.06]x[0.48,0.52]")));
  cout << "izls.size()=" << izls.size() << endl;
  
  Geometry::Rectangle<R> nr;
  Geometry::Zonotope<R> nz;
  Geometry::Zonotope<I,R> nez;
  Geometry::Zonotope<I> niz;
  Geometry::ListSet< Zonotope<I,R> > nezls;
  Geometry::ListSet< Zonotope<I,I> > nizls;
  
  Float x0=0;
  Float x1=0.4;
  Interval<Float> ivl1(0.4);
  Interval<Float> ivl0;
  
  time_type h(0.0625);
  cout << "h=" << h << endl;
  time_type t(0.25);
  cout << "t=" << t << endl;
  cout << endl;

  //Function evaluation sanity check
  cout << "vdp.image(" << r << ") = " << vdp.image(r) << endl;
  cout << "vdp.jacobian(" << r << ") = " << vdp.jacobian(r) << endl;
  cout << endl;
  
  // Integration step
  //nr=lohner.integration_step(vdp,r,h);
  //cout << nr << endl;
  LohnerIntegrator<R> plugin;
  VectorFieldEvolver<R> evolver(parameters,plugin);
  nez=evolver.integration_step(vdp,ez,h);
  cout << nez << endl << endl;
  cout << endl << endl;
  

  
  nezls=evolver.lower_integrate(vdp,ezls,t);
  cout << nezls << endl << endl;
  
  nezls=evolver.lower_reach(vdp,ezls,t);
  cout << nezls << endl << endl;
  
  // Affine vector field
  VectorFieldInterface<R>& avfr=avf;
  //AffineVectorField<R>& avfr=avf;
  nez=evolver.integration_step(avfr,ez,h);
  cout << nz << endl;
  cout << endl;
  
  return 0;
}


template<class R> 
int 
test_vector_field_evolver()
{
  typedef Interval<R> I;
  cout << __PRETTY_FUNCTION__ << endl;
  
  EvolutionParameters<R> parameters;
  parameters.set_maximum_basic_set_radius(0.125);
  parameters.set_lock_to_grid_time(1.0);
  parameters.set_maximum_step_size(0.25);
  
  AffineIntegrator<R> plugin;
  VectorFieldEvolver<R> evolver(parameters,plugin);

  AffineVectorField<R> affine_vector_field(Matrix<R>("[-2,-1;1,-2]"),Vector<R>("[0.125,0.25]"));
  
  Rectangle<R> bb("[-4,4]x[-4,4]");
  Rectangle<R> r("[-3.125,-2.875]x[-0.125,0.125]");
  FiniteGrid<R> fg(bb,128);
  const Grid<R>& g(fg.grid());
  
  GridMaskSet<R> initial_set(fg);
  initial_set.adjoin(over_approximation(r,g));
  GridMaskSet<R> bounding_set(fg);
  bounding_set.adjoin(over_approximation(bb,g));
  
  //GridMaskSet<R> chainreach=affine.chainreach(affine_vector_field,initial_set,bounding_set);
  time_type integration_time=3;
  uint n=12;
  GridMaskSet<R> integrate_set=initial_set;
  GridMaskSet<R> found_set=initial_set;
  for(uint i=0; i!=n; ++i) {
    found_set=evolver.integrate(affine_vector_field,found_set,bounding_set,time_type(integration_time/n));
    integrate_set.adjoin(found_set);
  }
  GridMaskSet<R> reach_set=evolver.reach(affine_vector_field,initial_set,bounding_set,integration_time);
  cout << endl;
  
  // set_integrator_verbosity(4);
  
  integration_time=0.5;
  evolver.parameters().set_grid_length(0.0625);
  PolyhedralSet<R> polyhedral_initial_set=PolyhedralSet<R>(Matrix<R>("[-2,0;0,-1;1,1]"),Vector<R>("[-1,-1,3]"));
  SetInterface<R>* polyhedral_initial_set_ptr=&polyhedral_initial_set;
  SetInterface<R>* polyhedral_integrate_set_ptr=evolver.integrate(affine_vector_field,*polyhedral_initial_set_ptr,integration_time);
  SetInterface<R>* polyhedral_reach_set_ptr=evolver.reach(affine_vector_field,*polyhedral_initial_set_ptr,integration_time);

  cout << "polyhedral_initial_set=" << *polyhedral_initial_set_ptr << endl;
  cout << "polyhedral_integrate_set=" << *polyhedral_integrate_set_ptr << endl;
  cout << "polyhedral_reach_set=" << *polyhedral_reach_set_ptr << endl;

  cout << endl;

  //Grid<R> grid(Vector<R>("[0.125,0.125]"));
  Grid<R> grid(Vector<R>("[0.0625,0.0625]"));
  ListSet< Rectangle<R> > rectangle_list_initial_set=lower_approximation(*polyhedral_initial_set_ptr,grid);
  ListSet< Rectangle<R> > rectangle_list_integrate_set=evolver.integrate(affine_vector_field,rectangle_list_initial_set,integration_time);
  ListSet< Rectangle<R> > rectangle_list_reach_set=evolver.reach(affine_vector_field,rectangle_list_initial_set,integration_time);

  cout << rectangle_list_initial_set << endl;
  cout << rectangle_list_integrate_set << endl;
  cout << rectangle_list_reach_set << endl;

  ListSet< Zonotope<I,R> > zonotope_list_initial_set=rectangle_list_initial_set;
  ListSet< Zonotope<I,R> > zonotope_list_integrate_set=evolver.lower_integrate(affine_vector_field,zonotope_list_initial_set,integration_time);
  ListSet< Zonotope<I,R> > zonotope_list_reach_set=evolver.lower_reach(affine_vector_field,zonotope_list_initial_set,integration_time);

  cout << zonotope_list_initial_set << endl;
  cout << zonotope_list_integrate_set << endl;
  cout << zonotope_list_reach_set << endl;

  cout << endl;

  epsfstream eps;
  eps.open("test_vector_field_evolver-1.eps",bb);
  eps << fill_colour(green) << reach_set;
  eps << fill_colour(yellow) << integrate_set;
  eps << fill_colour(blue) << initial_set;
  eps.close();

  eps.open("test_vector_field_evolver-2.eps",bb);
  eps << line_style(false);
  eps << fill_colour(green) << dynamic_cast<ListSet< Zonotope<I,I> >&>(*polyhedral_reach_set_ptr);
  eps << fill_colour(yellow) << dynamic_cast<ListSet< Zonotope<I,I> >&>(*polyhedral_integrate_set_ptr);
  eps << fill_colour(blue) << *polyhedral_initial_set_ptr;
  eps.close();
  
  eps.open("test_vector_field_evolver-3.eps",bb);
  eps << line_style(false);
  eps << fill_colour(green) << rectangle_list_reach_set;
  eps << fill_colour(yellow) << rectangle_list_integrate_set;
  eps << fill_colour(blue) << rectangle_list_initial_set;
  eps << fill_colour(red) << *polyhedral_initial_set_ptr;
  eps.close();
  
  eps.open("test_vector_field_evolver-4.eps",bb);
  eps << line_style(false);
  eps << fill_colour(green) << zonotope_list_reach_set;
  eps << fill_colour(yellow) << zonotope_list_integrate_set;
  eps << fill_colour(blue) << zonotope_list_initial_set;
  eps.close();
  

  return 0;
}
