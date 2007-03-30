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
#include "evaluation/lohner_integrator.h"
#include "evaluation/affine_integrator.h"
#include "models/vanderpol.h"
#include "output/epsfstream.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;
using namespace Ariadne::Output;
using namespace std;

template<class R> int test_integration_step();
template<class R> int test_integrate();

int main() {
  test_integration_step<Float>();
  test_integrate<Float>();
  return 0;
}

template<class R> 
int 
test_integration_step()
{
  cout << __PRETTY_FUNCTION__ << endl;
  typedef Interval<R> I;

  // Test constructor/destructor
  Integrator<R>* integrator_ptr;
  LohnerIntegrator<R>* lohner_ptr;
  lohner_ptr=new LohnerIntegrator<R>(0.125,0.5,0.125);
  delete lohner_ptr;
  integrator_ptr=new LohnerIntegrator<R>(0.125,0.5,0.125);
  delete integrator_ptr;
  
  LohnerIntegrator<R> lohner=LohnerIntegrator<R>(0.125,0.5,0.125);
  
  AffineVectorField<R> avf=AffineVectorField<R>(Matrix<R>("[-0.25,-1.0;+1.0,-0.25]"),Vector<R>("[0.25,0.0]"));
  cout << "avf=" << avf << endl;
  VanDerPolEquation<R> vdp=VanDerPolEquation<R>(R(0.865));
  cout << "vpd=" << avf << endl;

  Rectangle<R> r=Rectangle<R>("[0.98,1.02]x[0.48,0.52]");
  cout << "r=" << r << endl;
  Zonotope<R> z=Zonotope<R>(r);
  cout << "z=" << z << endl;
  Zonotope<I> iz=Zonotope<R>(r);
  cout << "iz=" << iz << endl;

  ListSet< Zonotope<R> > zls=ListSet< Zonotope<R> >(z);
  zls.adjoin(Zonotope<R>(Rectangle<R>("[1.02,1.06]x[0.48,0.52]")));
  cout << "zls.size()=" << zls.size() << endl;
  
  ListSet< Zonotope<I> > izls=ListSet< Zonotope<I> >(iz);
  izls.adjoin(Zonotope<I>(Rectangle<R>("[1.02,1.06]x[0.48,0.52]")));
  cout << "izls.size()=" << izls.size() << endl;
  
  Geometry::Rectangle<R> nr;
  Geometry::Zonotope<R> nz;
  Geometry::Zonotope<I> niz;
  Geometry::ListSet< Zonotope<R> > nzls;
  Geometry::ListSet< Zonotope<I> > nizls;
  
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
  niz=lohner.integration_step(vdp,iz,h);
  cout << niz << endl << endl;
  cout << endl << endl;
  

  
  // Integrate
  //nr=lohner.integrate(vdp,r,t);
  //cout << nr << endl;
  niz=lohner.integrate(vdp,iz,t);
  cout << niz << endl << endl;;
  
  nizls=lohner.integrate(vdp,izls,t);
  cout << nizls << endl << endl;
  
  nizls=lohner.reach(vdp,izls,t);
  cout << nizls << endl << endl;
  
  // Affine vector field
  VectorField<R>& avfr=avf;
  //AffineVectorField<R>& avfr=avf;
  niz=lohner.integration_step(avfr,iz,h);
  cout << nz << endl;
  cout << endl;
  
  return 0;
}


template<class R> 
int 
test_integrate()
{
  cout << __PRETTY_FUNCTION__ << endl;
  
  AffineIntegrator<R> affine(0.0625,1.0,0.25);
  
  AffineVectorField<R> affine_vector_field(Matrix<R>("[-2,-1;1,-2]"),Vector<R>("[1,1]"));
  
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
    found_set=affine.integrate(affine_vector_field,found_set,bounding_set,time_type(integration_time/n));
    integrate_set.adjoin(found_set);
  }
  GridMaskSet<R> reach_set=affine.reach(affine_vector_field,initial_set,bounding_set,integration_time);
  
  epsfstream eps;
  eps.open("test_integrate.eps",bb);
  eps << reach_set;
  eps.set_fill_colour("blue");
  eps << integrate_set;
  eps.set_fill_colour("yellow");
  eps << initial_set;
  eps.close();
  cout << endl;

  return 0;
}
