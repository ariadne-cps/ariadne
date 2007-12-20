/***************************************************************************
 *            test_vector_field_evolver.cc
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
#include "geometry/empty_set.h"
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

int main(int nargs, const char* args[]) 
{
  int verbosity = 0;
  if(nargs>1) {
    verbosity = std::atoi(args[1]);
  }
  set_integrator_verbosity(verbosity);
  test_integrator<Flt>();
  test_vector_field_evolver<Flt>();
  cout << "INCOMPLETE " << flush;
  return 0;
}

template<class R> 
int 
test_integrator()
{
  cout << __PRETTY_FUNCTION__ << endl;
  typedef Interval<R> I;

  EvolutionParameters<R> parameters;
  parameters.set_maximum_basic_set_radius(0.25);
  parameters.set_maximum_step_size(0.125);
  
  // Test constructor/destructor
  LohnerIntegrator<R>* lohner_ptr;
  lohner_ptr=new LohnerIntegrator<R>();
  delete lohner_ptr;
  
  AffineVectorField<R> avf=AffineVectorField<R>(Matrix<R>("[-0.25,-1.0;+1.0,-0.25]"),Vector<R>("[0.25,0.0]"));
  cout << "avf=" << avf << endl;
  R mu=0.865;
  VanDerPolEquation<R> vdp=VanDerPolEquation<R>(Point<R>(1,&mu));
  cout << "vpd=" << avf << endl;

  Rectangle<R> r=Rectangle<R>("[0.98,1.02]x[0.48,0.52]");
  cout << "r=" << r << endl;
  Zonotope<R,UniformErrorTag> ez=Zonotope<R,UniformErrorTag>(r);
  cout << "ez=" << ez << endl;

  ListSet< Zonotope<R,UniformErrorTag> > ezls=ListSet< Zonotope<R,UniformErrorTag> >(ez);
  ezls.adjoin(Zonotope<R>(Rectangle<R>("[1.02,1.06]x[0.48,0.52]")));
  cout << "ezls.size()=" << ezls.size() << endl;
  
  Geometry::Rectangle<R> nr;
  Geometry::Zonotope<R> nz;
  Geometry::Zonotope<R,UniformErrorTag> nez;
  Geometry::ListSet< Zonotope<R,UniformErrorTag> > nezls;
  
  Flt x0=0;
  Flt x1=0.4;
  Interval<Flt> ivl1(0.4);
  Interval<Flt> ivl0;
  
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
  parameters.set_maximum_basic_set_radius(0.25);
  parameters.set_lock_to_grid_time(0.5);
  parameters.set_maximum_step_size(0.125);
  parameters.set_grid_length(0.125);
  
  AffineIntegrator<R> affine_integrator;
  LohnerIntegrator<R> lohner_integrator;
  VectorFieldEvolver<R> evolver(parameters,lohner_integrator);

  AffineVectorField<R> affine_vector_field(Matrix<R>("[-2,-1;1,-2]"),Vector<R>("[0.125,0.25]"));
  
  Rectangle<R> bb("[-4,4]x[-4,4]");
  Rectangle<R> r("[-3.125,-2.875]x[-0.125,0.125]");
  FiniteGrid<R> fg(bb,128);
  //FiniteGrid<R> fg(bb,64);
  //FiniteGrid<R> fg(bb,32);
  const Grid<R>& g(fg.grid());
  Grid<R> fine_grid(Vector<R>("[0.0625,0.0625]"));
 
  time_type integration_time=2;
  uint n=12;

  evolver.parameters().set_grid_length(0.0625);



  //PolyhedralSet<R> polyhedral_initial_set=PolyhedralSet<R>(Matrix<R>("[-2,0;0,-1;1,1]"),Vector<R>("[-1,-1,3]"));
  PolyhedralSet<R> polyhedral_initial_set=PolyhedralSet<R>(r);
  RectangularSet<R> rectangular_bounding_set=RectangularSet<R>(bb);
  SetInterface<R>* abstract_bounding_set_ptr=&rectangular_bounding_set;
  SetInterface<R>* abstract_initial_set_ptr=&polyhedral_initial_set;
  cout << "abstract_initial_set=" << *abstract_initial_set_ptr << endl;
  SetInterface<R>* abstract_integrate_set_ptr=evolver.integrate(affine_vector_field,*abstract_initial_set_ptr,time_type(integration_time/n));
  cout << "abstract_integrate_set=" << *abstract_integrate_set_ptr << endl;
  SetInterface<R>* abstract_final_set_ptr=evolver.integrate(affine_vector_field,*abstract_initial_set_ptr,integration_time);
  SetInterface<R>* abstract_reach_set_ptr=evolver.reach(affine_vector_field,*abstract_initial_set_ptr,integration_time);
  cout << "abstract_reach_set=" << *abstract_reach_set_ptr << endl;
  SetInterface<R>* abstract_chainreach_set_ptr=evolver.chainreach(affine_vector_field,*abstract_initial_set_ptr,*abstract_bounding_set_ptr);
  cout << "abstract_chainreach_set=" << *abstract_reach_set_ptr << endl;
  cout << endl;



  GridMaskSet<R> grid_initial_set(fg);
  grid_initial_set.adjoin(outer_approximation(polyhedral_initial_set,g));
  cout << "grid_initial_set.size()=" << grid_initial_set.size() << endl;
  GridMaskSet<R> grid_bounding_set(fg);
  grid_bounding_set.adjoin(over_approximation(bb,g));
  cout << "grid_bounding_set.size()=" << grid_bounding_set.size() << endl;
  GridMaskSet<R> grid_integrate_set=grid_initial_set;
  GridMaskSet<R> grid_found_set=grid_initial_set;
  for(uint i=0; i!=n; ++i) {
    grid_found_set=evolver.bounded_integrate(affine_vector_field,grid_found_set,grid_bounding_set,time_type(integration_time/n));
    grid_integrate_set.adjoin(grid_found_set);
  }
  cout << "grid_integrate_set.size()=" << grid_integrate_set.size() << endl;
  GridMaskSet<R> grid_final_set=grid_initial_set;
  grid_final_set=evolver.bounded_integrate(affine_vector_field,grid_initial_set,grid_bounding_set,integration_time/2);
  cout << "grid_final_set.size()=" << grid_final_set.size() << endl;
  cerr << "WARNING: VectorFieldEvolver::bounded_integrate(VectorFieldInterface,GridMaskSet,GridMaskSet,Rational) may have accuracy problems.\n";
  GridMaskSet<R> grid_reach_set=evolver.bounded_reach(affine_vector_field,grid_initial_set,grid_bounding_set,integration_time);
  cout << "grid_reach_set.size()=" << grid_reach_set.size() << endl;
  GridMaskSet<R> grid_chainreach_set=evolver.chainreach(affine_vector_field,grid_initial_set,grid_bounding_set);
  cout << "grid_chainreach_set.size()=" << grid_chainreach_set.size() << endl;



  //Grid<R> grid(Vector<R>("[0.125,0.125]"));
  ListSet< Rectangle<R> > rectangle_list_initial_set=point_approximation(polyhedral_initial_set,fine_grid);
  cout << "rectangle_list_initial_set.size()=" << rectangle_list_initial_set.size() << endl;
  ListSet< Rectangle<R> > rectangle_list_integrate_set=rectangle_list_initial_set;
  cout << "rectangle_list_integrate_set.size()=" << rectangle_list_integrate_set << endl;
  ListSet< Rectangle<R> > rectangle_list_found_set=rectangle_list_initial_set;
  for(uint i=0; i!=n; ++i) {
    rectangle_list_found_set=evolver.lower_integrate(affine_vector_field,rectangle_list_found_set,time_type(integration_time/n));
    rectangle_list_integrate_set.adjoin(rectangle_list_found_set);
  }
  ListSet< Rectangle<R> > rectangle_list_final_set=evolver.lower_integrate(affine_vector_field,rectangle_list_initial_set,integration_time);
  cout << "rectangle_list_final_set.size()=" << rectangle_list_final_set.size() << endl;
  ListSet< Rectangle<R> > rectangle_list_reach_set=evolver.lower_reach(affine_vector_field,rectangle_list_initial_set,integration_time);
  cout << "rectangle_list_reach_set.size()=" << rectangle_list_reach_set.size() << endl;

  cout << rectangle_list_initial_set << endl;
  cout << rectangle_list_final_set << endl;
  cout << rectangle_list_integrate_set << endl;
  cout << rectangle_list_reach_set << endl;
  cout << endl;

  ListSet< Zonotope<R,UniformErrorTag> > zonotope_list_initial_set=point_approximation(polyhedral_initial_set,fine_grid);
  ListSet< Zonotope<R,UniformErrorTag> > zonotope_list_integrate_set=zonotope_list_initial_set;
  ListSet< Zonotope<R,UniformErrorTag> > zonotope_list_found_set=zonotope_list_initial_set;
  for(uint i=0; i!=n; ++i) {
    zonotope_list_found_set=evolver.lower_integrate(affine_vector_field,zonotope_list_found_set,time_type(integration_time/n));
    zonotope_list_integrate_set.adjoin(zonotope_list_found_set);
  }
  ListSet< Zonotope<R,UniformErrorTag> > zonotope_list_final_set=evolver.lower_integrate(affine_vector_field,zonotope_list_initial_set,integration_time);
  ListSet< Zonotope<R,UniformErrorTag> > zonotope_list_reach_set=evolver.lower_reach(affine_vector_field,zonotope_list_initial_set,integration_time);
   
  cout << zonotope_list_initial_set << endl;
  cout << zonotope_list_integrate_set << endl;
  cout << zonotope_list_reach_set << endl;
  cout << endl;

  epsfstream eps;

  eps.open("test_vector_field_evolver-abstract.eps",bb);
  eps << line_style(false);
  eps << fill_colour(red) << *abstract_chainreach_set_ptr;
  eps << fill_colour(green) << *abstract_reach_set_ptr;
  eps << line_style(true);
  eps << fill_colour(cyan) << *abstract_integrate_set_ptr;
  eps << fill_colour(yellow) << *abstract_final_set_ptr;
  eps << fill_colour(blue) << *abstract_initial_set_ptr;
  eps << fill_colour(transparant) << polyhedral_initial_set;
  eps.close();
  
  eps.open("test_vector_field_evolver-grid.eps",bb);
  eps << fill_colour(red) << grid_chainreach_set;
  eps << fill_colour(green) << grid_reach_set;
  eps << fill_colour(cyan) << grid_integrate_set;
  eps << fill_colour(yellow) << grid_final_set;
  eps << fill_colour(blue) << grid_initial_set;
  eps.close();

  eps.open("test_vector_field_evolver-rlist.eps",bb);
  eps << line_style(false);
  eps << fill_colour(green) << rectangle_list_reach_set;
  eps << line_style(true);
  eps << fill_colour(cyan) << rectangle_list_integrate_set;
  eps << fill_colour(yellow) << rectangle_list_final_set;
  eps << fill_colour(blue) << rectangle_list_initial_set;
  eps << fill_colour(transparant) << polyhedral_initial_set;
  eps.close();
  
  eps.open("test_vector_field_evolver-zlist.eps",bb);
  eps << line_style(false);
  eps << fill_colour(green) << zonotope_list_reach_set;
  eps << line_style(true);
  eps << fill_colour(cyan) << zonotope_list_integrate_set;
  eps << fill_colour(yellow) << zonotope_list_final_set;
  eps << fill_colour(blue) << zonotope_list_initial_set;
  eps << fill_colour(transparant) << polyhedral_initial_set;
  eps.close();
  

  return 0;
}
