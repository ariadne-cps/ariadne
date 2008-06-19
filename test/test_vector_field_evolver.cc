/***************************************************************************
 *            test_vector_field_evolver.cc
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
#include "evaluation/standard_integrator.h"
#include "evaluation/standard_subdivider.h"
#include "evaluation/cascade_reducer.h"
#include "models/vanderpol.h"
#include "output/epsstream.h"
#include "output/logging.h"

#include "test.h"

using namespace std;
using namespace Ariadne;
using Ariadne::Models::VanDerPolEquation;

template<class R> 
class TestVectorFieldEvolver
{
 public:
  void test() const;
};

int main() 
{
  TestVectorFieldEvolver<Flt>().test();
  return ARIADNE_TEST_FAILURES;
}

template<class R> 
void TestVectorFieldEvolver<R>::test() const
{
  cout << __PRETTY_FUNCTION__ << endl;
  typedef Interval<R> I;

  EvolutionParameters<R> parameters;
  parameters.set_maximum_basic_set_radius(0.25);
  parameters.set_maximum_step_size(0.125);
  
  // Test constructor/destructor
  StandardIntegrator< Zonotope<R> > integrator;
  StandardSubdivider< Zonotope<R> > subdivider;
  CascadeReducer< Zonotope<R> > reducer(3);

  VectorFieldEvolver< Zonotope<R> > evolver(parameters,integrator,subdivider,reducer);
  
  AffineVectorField<R> avf=AffineVectorField<R>(Matrix<R>("[-0.25,-1.0;+1.0,-0.25]"),Vector<R>("[0.25,0.0]"));
  cout << "avf=" << avf << endl;
  R mu=0.865;
  VanDerPolEquation<R> vdp=VanDerPolEquation<R>(Point<R>(1,&mu));
  cout << "vpd=" << avf << endl;

  Box<R> r=Box<R>("[0.98,1.02]x[0.48,0.52]");
  cout << "r=" << r << endl;
  Zonotope<R> z(r);
  cout << "z=" << z << endl;

  Box<R> bounding_box=Box<R>("[-4,4]x[-4,4]") ;
  Box<R> eps_bounding_box=bounding_box.neighbourhood(0.1);
 
  time_type h(0.0625);
  cout << "h=" << h << endl;
  time_type t(0.25);
  cout << "t=" << t << endl;
  cout << endl;

  //Function evaluation sanity check
  cout << "vdp.evaluate(" << r << ") = " << vdp.evaluate(r) << endl;
  cout << "vdp.jacobian(" << r << ") = " << vdp.jacobian(r) << endl;
  cout << endl;
  
  Zonotope<R> initial_set = z;
  ListSet< Zonotope<R> > evolve_set = evolver.evolve(vdp,z,1);
  ListSet< Zonotope<R> > reach_set = evolver.reach(vdp,z,1);
  
  epsfstream eps;
  eps.open("test_vector_field_evolver.eps",eps_bounding_box);
  eps << line_style(true);
  eps << fill_colour(cyan) << reach_set;
  eps << fill_colour(yellow) << evolve_set;
  eps << fill_colour(blue) << initial_set;
  eps.close();
}
