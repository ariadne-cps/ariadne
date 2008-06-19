/***************************************************************************
 *            test_applicator.cc
 *
 *  Copyright  2005-8  Pieter Collins
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

#include "test_float.h"

#include "base/pointer.h"
#include "geometry/point.h"
#include "geometry/zonotope.h"
#include "geometry/polytope.h"
#include "evaluation/evolution_parameters.h"
#include "evaluation/standard_applicator.h"
#include "output/epsstream.h"
#include "output/logging.h"

#include "models/henon.h"

#include "test.h"

using namespace std;
using namespace Ariadne;
using Ariadne::Models::HenonMap;
using Ariadne::Models::HenonInverseMap;

template<class R> 
class TestApplicator
{
 public:
  void test() const;
};

template<class R> 
void TestApplicator<R>::test() const
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
  

  //Test evaluation on different classes of sets
  StandardApplicator<ZBS> zonotope_applicator;
  Zonotope<R> fz=zonotope_applicator.apply(henon,z);

  //StandardApplicator<PBS> polytope_applicator;
  //Polytope<R> fpl=polytope_applicator.apply(henon,pl);

  Zonotope<R> pfz=zonotope_applicator.apply(henon_inverse,fz);
  cout << "z=" << z << " fz=" << fz << " pfz="<< pfz << endl;
  
  epsfstream eps;

}

int main() {
  TestApplicator<Flt>().test();
  return ARIADNE_TEST_FAILURES;
  
}
