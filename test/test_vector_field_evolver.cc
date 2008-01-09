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
#include "base/tuple.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "geometry/zonotope.h"
#include "geometry/list_set.h"
#include "geometry/empty_set.h"
#include "system/affine_vector_field.h"
#include "evaluation/evolution_parameters.h"
#include "evaluation/vector_field_evolver.h"
#include "evaluation/lohner_integrator.h"
#include "evaluation/affine_integrator.h"
#include "evaluation/standard_approximator.h"
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

  Box<R> r=Box<R>("[0.98,1.02]x[0.48,0.52]");
  cout << "r=" << r << endl;
  Zonotope<R> z(r);
  cout << "z=" << z << endl;
  ConstraintSet<R> initial_set(r);
  cout << "initial_set=" << initial_set << endl;

 
  Geometry::Box<R> nr;
  
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
  cout << "vdp.evaluate(" << r << ") = " << vdp.evaluate(r) << endl;
  cout << "vdp.jacobian(" << r << ") = " << vdp.jacobian(r) << endl;
  cout << endl;
  
  // Integration step
  //nr=lohner.integration_step(vdp,r,h);
  //cout << nr << endl;
  StandardApproximator< Zonotope<R> > approximator;
  StandardBounder<R> bounder;
  LohnerIntegrator<R> integrator;
  AffineIntegrator<R> affine_integrator;
  VectorFieldEvolver< Zonotope<R> > evolver(parameters,integrator,approximator);
  Box<R> bb;
  make_lpair(h,bb)=bounder.flow_bounds(vdp,z.bounding_box(),h);
  Zonotope<R> nz=integrator.integration_step(vdp,z,h,bb);
  cout << nz << endl << endl;
  cout << endl << endl;
  

  
  SetInterface<R>* evolve_ptr=evolver.lower_evolve(vdp,initial_set,t);
  cout << *evolve_ptr << endl << endl;
  
  SetInterface<R>* reach_ptr=evolver.lower_reach(vdp,initial_set,t);
  cout << *reach_ptr << endl << endl;
  
  // Affine vector field
  VectorField<R>& avfr=avf;
  //AffineVectorField<R>& avfr=avf;
  nz=affine_integrator.integration_step(avfr,z,h,bb);
  cout << nz << endl;
  cout << endl;
  
  return 0;
}


template<class R> 
int 
test_vector_field_evolver()
{

  typedef Interval<R> I;
  typedef Zonotope<R> BS;
  cout << __PRETTY_FUNCTION__ << endl;
  
  EvolutionParameters<R> parameters;
  parameters.set_maximum_basic_set_radius(0.25);
  parameters.set_lock_to_grid_time(0.5);
  parameters.set_maximum_step_size(0.125);
  parameters.set_grid_length(0.125);
  
  AffineIntegrator<R> affine_integrator;
  LohnerIntegrator<R> lohner_integrator;
  VectorFieldEvolver<BS> evolver(parameters,lohner_integrator);

  AffineVectorField<R> affine_vector_field(Matrix<R>("[-2,-1;1,-2]"),Vector<R>("[0.125,0.25]"));
  
  Box<R> bb("[-4,4]x[-4,4]");
  Box<R> r("[-3.125,-2.875]x[-0.125,0.125]");
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
  Box<R> bounding_box=bb;
  SetInterface<R>* abstract_initial_set_ptr=&polyhedral_initial_set;
  cout << "abstract_initial_set=" << *abstract_initial_set_ptr << endl;
  SetInterface<R>* abstract_integrate_set_ptr=evolver.upper_evolve(affine_vector_field,*abstract_initial_set_ptr,time_type(integration_time/n));
  cout << "abstract_integrate_set=" << *abstract_integrate_set_ptr << endl;
  SetInterface<R>* abstract_final_set_ptr=evolver.upper_evolve(affine_vector_field,*abstract_initial_set_ptr,integration_time);
  SetInterface<R>* abstract_reach_set_ptr=evolver.upper_reach(affine_vector_field,*abstract_initial_set_ptr,integration_time);
  cout << "abstract_reach_set=" << *abstract_reach_set_ptr << endl;
  SetInterface<R>* abstract_chainreach_set_ptr=evolver.chainreach(affine_vector_field,*abstract_initial_set_ptr);
  cout << "abstract_chainreach_set=" << *abstract_reach_set_ptr << endl;
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
  

  return 0;
}
