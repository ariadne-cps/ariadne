/***************************************************************************
 *            test_integration_step.cc
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

#include "real_typedef.h"

#include "ariadne.h"
#include "numeric/numerical_types.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"
#include "geometry/list_set.h"
#include "system/affine_vector_field.h"
#include "evaluation/lohner_integrator.h"
#include "models/vanderpol.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;
using namespace std;

template<class R> int test_integration_step();

int main() {
  test_integration_step<Real>();
  return 0;
}

template<class R> 
int 
test_integration_step()
{
  
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
  Parallelotope<R> p=Parallelotope<R>(r);
  cout << "p=" << p << endl;
  Parallelotope< Interval<R> > ip=Parallelotope<R>(r);
  cout << "ip=" << ip << endl;
  Zonotope<R> z=Zonotope<R>(r);
  cout << "z=" << z << endl;
  Zonotope< Interval<R> > iz=Zonotope<R>(r);
  cout << "iz=" << iz << endl;

  ListSet<R,Parallelotope> pls=ListSet<R,Parallelotope>(p);
  cout << "pls.size()=" << pls.size() << endl;
  ListSet<R,Zonotope> zls=ListSet<R,Zonotope>(z);
  zls.adjoin(Zonotope<R>(Rectangle<R>("[1.02,1.06]x[0.48,0.52]")));
  cout << "zls.size()=" << zls.size() << endl;
  
  Geometry::Rectangle<R> nr;
  Geometry::Parallelotope<R> np;
  Geometry::Parallelotope< Interval<R> > nip;
  Geometry::Zonotope<R> nz;
  Geometry::Zonotope< Interval<R> > niz;
  Geometry::ListSet<R,Parallelotope> npls;
  Geometry::ListSet<R,Zonotope> nzls;
  
  Real x0=0;
  Real x1=0.4;
  Interval<Real> ivl1(0.4);
  Interval<Real> ivl0;
  
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
  nz=lohner.integration_step(vdp,z,h);
  cout << nz << endl << endl;
  nip=lohner.integration_step(vdp,ip,h);
  cout << nip << endl << endl;
  np=lohner.integration_step(vdp,p,h);
  cout << np << endl << endl;
  cout << endl << endl;
  

  
  // Integrate
  //nr=lohner.integrate(vdp,r,t);
  //cout << nr << endl;
  nz=lohner.integrate(vdp,z,t);
  cout << nip << endl << endl;;
  np=lohner.integrate(vdp,p,t);
  cout << np << endl << endl;;
  
  // 'integrate' not defined for fuzzy sets
  //nip=lohner.integrate(vdp,ip,t);
  //niz=lohner.integrate(vdp,iz,t);

  //npls=lohner.integrate(vdp,pls,t);
  //cout << npls << endl << endl;;
  nzls=lohner.integrate(vdp,zls,t);
  cout << nzls << endl << endl;
  
  nzls=lohner.reach(vdp,zls,t);
  cout << nzls << endl << endl;
  
  // Affine vector field
  VectorField<R>& avfr=avf;
  //AffineVectorField<R>& avfr=avf;
  np=lohner.integration_step(avfr,p,h);
  cout << np << endl;
  nip=lohner.integration_step(avfr,ip,h);
  cout << nip << endl;
  nz=lohner.integration_step(avfr,z,h);
  cout << nz << endl;
  niz=lohner.integration_step(avfr,iz,h);
  cout << niz << endl;

  
  return 0;
}
