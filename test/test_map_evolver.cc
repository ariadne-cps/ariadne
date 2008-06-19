/***************************************************************************
 *            test_map_evolver.cc
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
#include "geometry/rectangular_set.h"
#include "system/grid_multimap.h"
#include "evaluation/standard_applicator.h"
#include "evaluation/standard_subdivider.h"
#include "evaluation/cascade_reducer.h"
#include "evaluation/map_evolver.h"
#include "output/epsstream.h"
#include "output/logging.h"

#include "models/henon.h"

#include "test.h"

using namespace std;
using namespace Ariadne;
using Ariadne::Models::HenonMap;
using Ariadne::Models::HenonInverseMap;


template<class R> 
class TestMapEvolver
{
 public:
  TestMapEvolver();
  void test() const;
};


template<class R> 
TestMapEvolver<R>::TestMapEvolver()
{
}

template<class R>
void TestMapEvolver<R>::test() const
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
  cout << henon << endl << henon_inverse << endl << endl;

  Box<R> bounding_box=Box<R>("[-4,4]x[-4,4]") ;
  Box<R> eps_bounding_box=bounding_box.neighbourhood(0.1);
  
  Box<R> bx=Box<R>("[1.49,1.51]x[0.49,0.51]");
  Zonotope<R> z(bx);
  Polytope<R> pl(bx);
  
  EvolutionParameters<R> parameters;
  StandardApplicator< Zonotope<R> > applicator;
  StandardSubdivider< Zonotope<R> > subdivider;
  CascadeReducer< Zonotope<R> > reducer(3);
  
  Evolver< Map<R>, Zonotope<R> > evolver(parameters,applicator,subdivider,reducer);
  
  //Test evaluation on different classes of sets
  uint steps=3;
  ListSet< Zonotope<R> > fzl=evolver.reach(static_cast<const Map<R>&>(henon),z,Integer(steps));
  ARIADNE_ASSERT(fzl.size()==(steps+1u));
  cout << "z=" << z << "\nfzl=" << fzl << endl;
  Zonotope<R> fz=fzl[steps];
  cout << "fz=" << fz << endl;

  ListSet< Zonotope<R> > pfzl=evolver.reach(static_cast<const Map<R>&>(henon_inverse),fz,Integer(steps));
  cout << "pfzl=" << pfzl << endl;
  ARIADNE_ASSERT(pfzl.size()==(steps+1u));
  Zonotope<R> pfz=pfzl[steps];
  cout << "z=" << z << " fz=" << fz << " pfz="<< pfz << endl;
  
  Zonotope<R> initial_set = z;
  cout << "initial_set=" << initial_set << endl;
  ListSet< Zonotope<R> > evolve_set(fz);
  cout << "evolve_set=" << evolve_set << endl;
  ListSet< Zonotope<R> > reach_set(pfz);
  cout << "reach_set=" << reach_set << endl;
  

  epsfstream eps;
  eps.open("test_map_evolver.eps",eps_bounding_box);
  eps << line_style(true);
  eps << fill_colour(cyan) << reach_set;
  eps << fill_colour(yellow) << evolve_set;
  eps << fill_colour(blue) << initial_set;
  eps.close();



}

int main() {
  TestMapEvolver<Flt>().test();
  return ARIADNE_TEST_FAILURES;
}

